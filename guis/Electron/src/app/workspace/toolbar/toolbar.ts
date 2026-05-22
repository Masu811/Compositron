import { Component } from "@angular/core";
import { Router } from "@angular/router";
import { Dtype, MultiCampaign, DtypeSelection, parseSelection, PlotTrace } from "../../types";
import { AppData } from "../../app_data";

@Component({
  selector: "toolbar",
  templateUrl: "toolbar.html",
  styleUrl: "toolbar.css",
})
export class Toolbar {
  constructor(private router: Router, public appData: AppData) { }

  async importData() {
    const path = await (window as any).api.showOpenDialog();

    if (path == null) {
      return;
    }

    const response = await fetch(`http://127.0.0.1:8000/import_data?path=${path}`);
    const newData = await response.json() as Array<MultiCampaign>;

    this.appData.data.set(newData);
    this.appData.projectLoaded.set(true);
  }

  async deleteData() {
    const payload = parseSelection(this.appData.selection());

    const response = await fetch("http://127.0.0.1:8000/delete_data", {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify(payload),
    });

    const newData = await response.json() as Array<MultiCampaign>;

    this.appData.deselect.update(val => !val);
    this.appData.clearSelection();
    this.appData.data.set(newData);
  }

  async fetchMetadata(kind: Dtype) {
    const idcs = parseSelection(this.appData.selection())[kind][0];

    const response = await fetch(`http://127.0.0.1:8000/getMetadata/${idcs}`);
    if (!response.ok) {
      console.error(`Failed to fetch metadata for indices ${idcs}`);
      return;
    }
    const data = await response.json();

    this.appData.metadata.set(data);
    this.router.navigate(["/metadata"]);
  }

  analyzeSingles() {
    this.appData.dialogType.set("singleAnalyze");
    this.appData.dialogOpen.set(true);
  }

  async showSpectrum(kind: Dtype) {
    const idcs = parseSelection(this.appData.selection())[kind][0];

    let endpoint;

    switch (kind) {
      case Dtype.S:
        endpoint = "showSingleSpectrum";
        break;
      case Dtype.C:
        endpoint = "showCoincidenceSpectrum";
        break;
      case Dtype.M:
        endpoint = "showSingles";
        break;
      default:
        console.error(`Unsupported Dtype ${kind} supplied to function showSpectrum`);
        return;
    }

    const response = await fetch(`http://127.0.0.1:8000/${endpoint}/${idcs}`);
    if (!response.ok) {
      console.error(`Failed to fetch metadata for indices ${idcs}`);
      return;
    }
    const data = await response.json();

    if (kind == Dtype.C) {
      //@ts-ignore
      data.map((trace: PlotTrace) => trace.z = trace.z.map(row => row.map(i => Math.log(i))));
    }

    this.appData.plotData.set(data);
    this.router.navigate(["/plots"]);
  }

  print() {
    this.appData.dialogType.set("print");
    this.appData.dialogOpen.set(true);
  }

  parse() {
    this.appData.dialogType.set("parse");
    this.appData.dialogOpen.set(true);
  }

  filter() {
    this.appData.dialogType.set("filter");
    this.appData.dialogOpen.set(true);
  }

  sort() {
    this.appData.dialogType.set("sort");
    this.appData.dialogOpen.set(true);
  }

  split() {
    this.appData.dialogType.set("split");
    this.appData.dialogOpen.set(true);
  }

  async merge() {
    const payload = parseSelection(this.appData.selection());

    const response = await fetch("http://127.0.0.1:8000/merge", {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify(payload),
    });

    if (!response.ok) {
      console.error(response);
      return;
    }

    const newData = await response.json() as Array<MultiCampaign>;

    this.appData.deselect.update(val => !val);
    this.appData.clearSelection();
    this.appData.data.set(newData);
  }

  cluster() {
    this.appData.dialogType.set("cluster");
    this.appData.dialogOpen.set(true);
  }

  average() {
    this.appData.dialogType.set("average");
    this.appData.dialogOpen.set(true);
  }

  sum() {
    this.appData.dialogType.set("sum");
    this.appData.dialogOpen.set(true);
  }

  plot() {
    this.appData.dialogType.set("plot");
    this.appData.dialogOpen.set(true);
  }

  ratio_plotter() {
    this.appData.dialogType.set("ratio_plotter");
    this.appData.dialogOpen.set(true);
  }
}
