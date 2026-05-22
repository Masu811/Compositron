from enum import Enum
from typing import Any, TypedDict, Callable, Literal
from dataclasses import dataclass, asdict

import numpy as np
from numpy.typing import NDArray
from fastapi import FastAPI, WebSocket, WebSocketDisconnect
from fastapi.responses import JSONResponse

from stacs.coinc import CoincidenceSpectrum
from stacs.single import SingleSpectrum
from stacs import MultiCampaign, MeasurementCampaign, DopplerMeasurement, measurement

from parsers import Parser, parse_parser


data: list[MultiCampaign] = []


app = FastAPI()


@app.websocket("/health")
async def health(websocket: WebSocket):
    await websocket.accept()
    try:
        while True:
            await websocket.receive()
    except (WebSocketDisconnect, RuntimeError):
        print("Client disconnected")
        pass


def tree(m: MultiCampaign | MeasurementCampaign | DopplerMeasurement):
    if isinstance(m, DopplerMeasurement):
        return {
            "name": m.name,
            "singles": [s.detname for s in m.singles],
            "coinc": [c.detpair for c in m.coinc],
        }

    if isinstance(m, MeasurementCampaign):
        return {
            "name": m.name,
            "measurements": [tree(n) for n in m],
        }

    return {
        "name": m.name,
        "campaigns": [tree(n) for n in m]
    }


@app.get("/import_data")
def import_data(path: str):
    mult = MultiCampaign(
        path,
        autocompute_singles=False,
        autocompute_coinc=False,
    )
    data.append(mult)
    return JSONResponse([tree(t) for t in data])


class Dtype(Enum):
    MULT = "MULT"
    MC = "MC"
    M = "M"
    S = "S"
    C = "C"


Selection = dict[Dtype, list[str]]


class ParsedSelection(TypedDict):
    MULT: list[tuple[int]]
    MC: list[tuple[int, int]]
    M: list[tuple[int, int, int]]
    S: list[tuple[int, int, int, str]]
    C: list[tuple[int, int, int, str]]


ParsedIndices = tuple[int] | tuple[int, int] | tuple[int, int, int] | tuple[int, int, int, str]


selection_types = [int, int, int, str]


def parse_single_selection(idcs: str) -> ParsedIndices:
    fields = idcs.split("-")
    parsed_idcs = []

    for field, type in zip(fields, selection_types):
        parsed_idcs.append(type(field))

    return tuple(parsed_idcs)


def parse_selection(selection: Selection) -> ParsedSelection:
    parsed_selection: ParsedSelection = {
        "MULT": [],
        "MC": [],
        "M": [],
        "S": [],
        "C": [],
    }

    for dtype, idx_arr in selection.items():
        for str_idcs in idx_arr:
            parsed_selection[dtype.value].append(parse_single_selection(str_idcs))

    return parsed_selection


@app.post("/delete_data")
def delete_data(selection: Selection):
    parsed_selection = parse_selection(selection)

    remove_mult: list[tuple[int]] = sorted(
        parsed_selection["MULT"],
        key=lambda x: x[0],
        reverse=True,
    )

    remove_mc: list[tuple[int, int]] = sorted(
        filter(
            lambda x: (
                not any(x[0] == y[0] for y in remove_mult)
            ),
           parsed_selection["MC"]
        ),
        key=lambda x: x[1],
        reverse=True,
    )

    remove_m: list[tuple[int, int, int]] = sorted(
        filter(
            lambda x: (
                not any(x[0] == y[0] for y in remove_mult)
                and not any(x[0] == y[0] and x[1] == y[1] for y in remove_mc)
            ),
           parsed_selection["M"]
        ),
        key=lambda x: x[2],
        reverse=True,
    )

    remove_s: list[tuple[int, int, int, str]] = sorted(
        filter(
            lambda x: (
                not any(x[0] == y[0] for y in remove_mult)
                and not any(x[0] == y[0] and x[1] == y[1] for y in remove_mc)
                and not any(x[0] == y[0] and x[1] == y[1] and x[2] == y[2] for y in remove_m)
            ),
           parsed_selection["S"]
        )
    )

    remove_c: list[tuple[int, int, int, str]] = sorted(
        filter(
            lambda x: (
                not any(x[0] == y[0] for y in remove_mult)
                and not any(x[0] == y[0] and x[1] == y[1] for y in remove_mc)
                and not any(x[0] == y[0] and x[1] == y[1] and x[2] == y[2] for y in remove_m)
            ),
           parsed_selection["C"]
        )
    )

    for idcs in remove_s:
        m = data[idcs[0]][idcs[1]][idcs[2]]
        s = m[idcs[3]]
        m.singles.remove(s)

    for idcs in remove_c:
        m = data[idcs[0]][idcs[1]][idcs[2]]
        c = m[idcs[3]]
        m.coinc.remove(c)

    for idcs in remove_m:
        data[idcs[0]][idcs[1]].measurements.pop(idcs[2])

    for idcs in remove_mc:
        data[idcs[0]].campaigns.pop(idcs[1])

    for idcs in remove_mult:
        data.pop(idcs[0])

    return JSONResponse([tree(t) for t in data])


def scalarize(x):
    return str(np.array(x).tolist())


def mult_to_json(x: MultiCampaign):
    vals = [
        "name", "path", "directory"
    ]
    return {"__name__": "MultiCampaign"} | {
        val: scalarize(getattr(x, val)) for val in vals
    } | {"length": len(x)}


def mc_to_json(x: MeasurementCampaign):
    vals = [
        "name", "path", "directory"
    ]
    return {"__name__": "MeasurementCampaign"} | {
        val: scalarize(getattr(x, val)) for val in vals
    } | {"length": len(x)}


def m_to_json(x: DopplerMeasurement):
    vals = [
        "name", "path", "directory", "filename", "filetype"
    ]
    return {"__name__": "DopplerMeasurement"} | {
        val: scalarize(getattr(x, val)) for val in vals
    } | {key: str(val) for key, val in x.metadata.items()}


def s_to_json(x: SingleSpectrum):
    vals = [
        "detname", "ecal", "eres", "peak_center", "peak_counts", "dpeak_counts",
        "counts", "dcounts", "s", "ds", "w", "dw", "v2p", "dv2p",
    ]
    return {"__name__": "SingleSpectrum"} | {
        val: scalarize(getattr(x, val)) for val in vals
    }


def c_to_json(x: CoincidenceSpectrum):
    vals = [
        "detpair", "ecal", "eres", "peak_counts", "dpeak_counts",
        "counts", "dcounts", "s", "ds", "w", "dw", "roi_params"
    ]
    return {"__name__": "CoincidenceSpectrum"} | {
        val: scalarize(getattr(x, val)) for val in vals
    }


def to_json(
    x: MultiCampaign
    | MeasurementCampaign
    | DopplerMeasurement
    | SingleSpectrum
    | CoincidenceSpectrum
):
    if isinstance(x, MultiCampaign):
        return mult_to_json(x)
    if isinstance(x, MeasurementCampaign):
        return mc_to_json(x)
    if isinstance(x, DopplerMeasurement):
        return m_to_json(x)
    if isinstance(x, SingleSpectrum):
        return s_to_json(x)
    return c_to_json(x)


@app.get("/getMetadata/{idcs}")
def get_metadata(idcs: str):
    parsed_idcs = parse_single_selection(idcs)

    x = data
    for i in parsed_idcs:
        x = x[i]

    return JSONResponse(to_json(x))


@dataclass
class SingleAnalyzeArgs:
    s_width: float
    w_width: float
    w_dist: float
    w_rightonly: bool
    peak_width: float
    bg_frac: float
    bg_corr: bool
    v2p_bounds: list[float]
    follow_peak_order: int


@app.post("/single_analyze")
def single_analyze(selection: Selection, args: SingleAnalyzeArgs):
    parsed_idcs = parse_selection(selection)

    for idcs in parsed_idcs["S"]:
        s = data[idcs[0]][idcs[1]][idcs[2]][idcs[3]]

        kwargs = asdict(args)

        order = kwargs.pop("follow_peak_order")

        kwargs["follow_peak"] = (
            "0th_order_correction" if order == 0 else
            "1st_order_correction" if order == 1 else
            ""
        )

        s.analyze(**kwargs)


@app.get("/showSingleSpectrum/{idcs}")
def show_SingleSpectrum(idcs: str):
    parsed_idcs = parse_single_selection(idcs)

    x = data
    for i in parsed_idcs:
        x = x[i]

    return JSONResponse([{
        "x": x.energies.tolist(),
        "y": x.spectrum.tolist(),
        "name": x.detname,
        "type": "scatter",
    }])


@app.get("/showCoincidenceSpectrum/{idcs}")
def show_CoincidenceSpectrum(idcs: str):
    parsed_idcs = parse_single_selection(idcs)

    x = data
    for i in parsed_idcs:
        x = x[i]

    return JSONResponse([{
        "x": (x.ecal[0][0] + x.ecal[0][1] * np.arange(x.window[0][0], x.window[0][1])).tolist(),
        "y": (x.ecal[1][0] + x.ecal[1][1] * np.arange(x.window[1][0], x.window[1][1])).tolist(),
        "z": x.hist.tolist(),
        "name": x.detpair,
        "type": "heatmap",
    }])


@app.get("/showSingles/{idcs}")
def show_singles(idcs: str):
    parsed_idcs = parse_single_selection(idcs)

    x = data
    for i in parsed_idcs:
        x = x[i]

    for s in x.singles:
        if s.energies is None:
            s.analyze()

    return JSONResponse([
        {
            "x": s.energies.tolist(),
            "y": s.spectrum.tolist(),
            "name": s.detname,
            "type": "scatter",
        }
        for s in x.singles
    ])


@dataclass
class Param:
    param: str
    target: str
    parser: Parser


class ParsedParams(TypedDict):
    params: list[str]
    targets: list[str | None]
    parsers: list[Callable]


def parse_targets_and_parsers(params: list[Param]) -> ParsedParams:
    return {
        "params": [p.param for p in params],
        "targets": [
            "doppler" if p.target.startswith("D")
            else "single" if p.target.startswith("S")
            else "coinc" if p.target.startswith("C")
            else None
            for p in params
        ],
        "parsers": [
            parse_parser(p.parser) for p in params
        ]
    }


@dataclass
class GenericArgs:
    params: list[Param]


@app.post("/print")
def print_data(selection: Selection, args: GenericArgs):
    kwargs = parse_targets_and_parsers(args.params)
    return [
        [
            {key: list(map(str, value)) for key, value in mc_data.items()}
            for mc_data in mult.get(**kwargs, use_pd=False)
        ] for mult in data
    ]


@app.post("/parse")
def parse(selection: Selection, args: GenericArgs):
    kwargs = parse_targets_and_parsers(args.params)
    for mult in data:
        mult.parse(**kwargs)


@dataclass
class FilterArgs:
    param: str
    target: str
    parser: Parser
    byValue: bool
    values: list[str]
    min: float
    max: float
    negative: bool


def filter_empty_mults():
    filtered_data = list(filter(lambda mult: len(mult) > 0, data))
    data.clear()
    data.extend(filtered_data)


def filter_empty_mcs():
    for mult in data:
        mult.campaigns = list(filter(lambda mc: len(mc) > 0, mult.campaigns))


def filter_empty_ms():
    for mult in data:
        for mc in mult:
            mc.measurements = list(filter(
                lambda m: len(m.singles) > 0 or len(m.coinc) > 0,
                mc.measurements
            ))


@app.post("/filter")
def filter_data(selection: Selection, args: FilterArgs):
    kwargs = {
        "values": args.values if args.byValue else None,
        "min": args.min,
        "max": args.max,
        "negative": args.negative,
    }

    kwargs.update(parse_targets_and_parsers([Param(
        param = args.param,
        target = args.target,
        parser = args.parser,
    )]))

    for mult in data:
        mult.filter(**kwargs)

    filter_empty_ms()
    filter_empty_mcs()
    filter_empty_mults()

    return JSONResponse([tree(t) for t in data])


@app.post("/sort")
def sort(selection: Selection, args: GenericArgs):
    kwargs = parse_targets_and_parsers(args.params)
    for mult in data:
        mult.sort(**kwargs)

    return JSONResponse([tree(t) for t in data])


@app.post("/split")
def split(selection: Selection, args: GenericArgs):
    kwargs = parse_targets_and_parsers(args.params)
    for mult in data:
        mult.split(**kwargs)

    return JSONResponse([tree(t) for t in data])


def safe_numpy(arr: NDArray):
    if isinstance(arr, np.generic):
        return None if np.isnan(arr) else arr.item()
    elif isinstance(arr, str):
        return arr
    elif hasattr(arr, "__iter__"):
        return [safe_numpy(elem) for elem in arr]
    else:
        return None if np.isnan(arr) else arr


@app.post("/get")
def get(selection: Selection, param: Param):
    kwargs = parse_targets_and_parsers([param])
    parsed_selection = parse_selection(selection)

    mult_selected = len(parsed_selection["MULT"]) > 0

    if mult_selected:
        selected_data = data[parsed_selection["MULT"][0][0]]
    else:
        selected_data = MultiCampaign(
            campaigns=[
                data[mc_idcs[0]][mc_idcs[1]]
                for mc_idcs in parsed_selection["MC"]
            ]
        )

    return JSONResponse([
        {key: safe_numpy(values) for key, values in foo.items()}
        for foo in selected_data.get(
            **kwargs,
            use_pd=False,
        )
    ])


@app.post("/merge")
def merge(selection: Selection):
    parsed_selection = parse_selection(selection)

    merged_mc = MultiCampaign(
        campaigns=[
            data[mc_idcs[0]].campaigns.pop(mc_idcs[1])
            for mc_idcs in sorted(parsed_selection["MC"], reverse=True)
        ]
    ).merge()

    if len(parsed_selection["MULT"]) == 1:
        data[parsed_selection["MULT"][0][0]].campaigns.append(merged_mc)
    elif all(mc_idcs[0] == parsed_selection["MC"][0][0] for mc_idcs in parsed_selection["MC"]):
        data[parsed_selection["MC"][0][0]].campaigns.append(merged_mc)
    else:
        data.append(MultiCampaign(campaigns=[merged_mc]))

    filter_empty_mults()

    return JSONResponse([tree(t) for t in data])


@dataclass
class ClusterArgs:
    params: list[Param]
    method: str
    clusters: list[list[float | None]]
    n_clusters: None | int


@app.post("/cluster")
def cluster(selection: Selection, args: ClusterArgs):
    parsed_selection = parse_selection(selection)
    kwargs = asdict(args)

    kwargs.pop("params")
    kwargs.update(parse_targets_and_parsers(args.params))

    clusters = kwargs.pop("clusters")
    n_clusters = kwargs.pop("n_clusters")
    method = kwargs.pop("method")

    if method == "n_clusters":
        kwargs["n_clusters"] = n_clusters
    elif method == "clusters":
        kwargs["clusters"] = np.array(clusters, ndmin=2).T
    elif method == "guess":
        kwargs["guess"] = np.array(clusters, ndmin=2).T

    for mc_idcs in parsed_selection["MC"]:
        data[mc_idcs[0]][mc_idcs[1]].cluster(**kwargs)


@dataclass
class AvgArgs:
    params: list[Param]
    keep_params: list[Param]


@app.post("/average")
def average(selection: Selection, args: AvgArgs):
    parsed_selection = parse_selection(selection)

    params = parse_targets_and_parsers(args.params)
    keep_params = parse_targets_and_parsers(args.keep_params)

    for mc_idcs in parsed_selection["MC"]:
        data[mc_idcs[0]][mc_idcs[1]].average(
            **params,
            keep_params = keep_params["params"],
            keep_targets = keep_params["targets"],
            keep_parsers = keep_params["parsers"],
        )

    return JSONResponse([tree(t) for t in data])


@app.post("/sum")
def sum(selection: Selection, args: AvgArgs):
    parsed_selection = parse_selection(selection)

    params = parse_targets_and_parsers(args.params)
    keep_params = parse_targets_and_parsers(args.keep_params)

    for mc_idcs in parsed_selection["MC"]:
        data[mc_idcs[0]][mc_idcs[1]].sum(
            **params,
            keep_params = keep_params["params"],
            keep_targets = keep_params["targets"],
            keep_parsers = keep_params["parsers"],
            analyze=False,
        )

    return JSONResponse([tree(t) for t in data])


@dataclass
class RatioPlotterArgs:
    reference: str | int
    comparands: None | list[str | int]
    normalize_eres: bool
    bins: float
    bg_sub: bool
    roi_width: None | float
    mirrored: bool


@app.post("/ratio_plotter")
def ratio_plotter(selection: Selection, args: RatioPlotterArgs):
    parsed_selection = parse_selection(selection)

    mc_idcs = parsed_selection["MC"][0]

    mc = data[mc_idcs[0]][mc_idcs[1]]

    energies, ratios, dratios = mc.calc_ratio(**asdict(args))

    return JSONResponse({
        "energies": {det: safe_numpy(arr) for det, arr in energies.items()},
        "ratios": {det: safe_numpy(arr) for det, arr in ratios.items()},
        "dratios": {det: safe_numpy(arr) for det, arr in dratios.items()},
    })
