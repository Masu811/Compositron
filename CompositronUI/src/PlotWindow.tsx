import { memo, useEffect, useState } from "react";
import { invoke } from "@tauri-apps/api/core";
import { open } from '@tauri-apps/plugin-dialog';
import uPlot from 'uplot';
import UplotReact from 'uplot-react';
import 'uplot/dist/uPlot.min.css';


async function importWaveforms(
  setFilepath: (filepath: string) => void,
  setNTraces: (nTraces: number) => void,
): Promise<void> {
  try {
    const path = await open({
      multiple: false,
      directory: false,
    });

    if (!path) return;

    setFilepath(path);

    setNTraces(await invoke("import_waveforms", { filepath: path }));
  } catch (e) {
    console.error(e);
  }
}


function ImporterButton({
    setFilepath, setNTraces
  }: {
    setFilepath: (filepath: string) => void,
    setNTraces: (nTraces: number) => void,
  }) {
  return (
    <div>
      <button onClick={() => importWaveforms(setFilepath, setNTraces)}>Import Data</button>
    </div>
  )
}


async function getTraceData(
  i: number,
  setTraceData: (traceData: Array<number>) => void,
  setPlotTitle: (plotTitle: string) => void,
) {
  try {
    setTraceData(await invoke("get_trace_data", { index: i }));
    setPlotTitle(`Trace ${i}`);
  } catch (e) {
    console.error(e);
  }
}


const TraceSelector = memo(({
  nTraces, setTraceData, setPlotTitle
}: {
    nTraces: number,
    setTraceData: (traceData: Array<number>) => void,
    setPlotTitle: (plotTitle: string) => void,
}) => {
  const traces = [...Array(nTraces).keys()].map((_, i) => {
    return <li key={i} onClick={() => getTraceData(i, setTraceData, setPlotTitle)}>{i}</li>;
  });
  return (
    <div className="trace-selector">
      <ul>
        {traces}
      </ul>
    </div>
  );
});


function PlotPanel({
  filepath, nTraces
}: {
  filepath: string,
  nTraces: number,
}) {
  const [traceData, setTraceData] = useState(Array<number>());
  const [plotTitle, setPlotTitle] = useState("Select a Trace");

  const [container, setContainer] = useState<HTMLDivElement | null>(null);
  const [size, setSize] = useState({ width: 1000, height: 800 });

  useEffect(() => {
    if (!container) return;

    const observer = new ResizeObserver(([entry]) => {
      setSize({
        width: entry.contentRect.width - 50,
        height: entry.contentRect.height - 50,
      });
    });

    observer.observe(container);

    return () => observer.disconnect();
  }, [container]);

  let opts = {
    title: plotTitle,
    id: "chart1",
    class: "my-chart",
    ...size,
    scales: { x: { time: false } },
    series: [
      {},
      {
        show: true,
        spanGaps: false,
        stroke: "steelblue",
        width: 1.5,
      }
    ],
  };

  if (!filepath) {
    return <h1>No Data Loaded</h1>;
  }

  return (
    <div id="plot-panel">
      <TraceSelector nTraces={nTraces} setTraceData={setTraceData} setPlotTitle={setPlotTitle}/>
      <div id="chart-container" ref={setContainer}>
        <UplotReact
          options={opts}
          data={[traceData.map((_, i) => i), traceData] as uPlot.AlignedData}
        />
      </div>
    </div>
  );
}


export function PlotWindow() {
  const [nTraces, setNTraces] = useState(0);
  const [filepath, setFilepath] = useState("");

  return (
    <div className="plot-window">
      <ImporterButton setFilepath={setFilepath} setNTraces={setNTraces} />
      <PlotPanel filepath={filepath} nTraces={nTraces}></PlotPanel>
    </div>
  );
}
