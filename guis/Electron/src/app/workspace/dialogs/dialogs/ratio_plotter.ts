import { Component, inject } from "@angular/core";
import { parseSelection, PlotAxes, PlotTrace } from "../../../types";
import { AppData } from "../../../app_data";
import { NonNullableFormBuilder, ReactiveFormsModule, Validators } from "@angular/forms";
import { Router } from "@angular/router";

@Component({
  templateUrl: "ratio_plotter.html",
  styleUrl: "../dialog.css",
  imports: [ReactiveFormsModule],
})
export class RatioPlotterDialog {
  private formBuilder = inject(NonNullableFormBuilder);

  formModel = this.formBuilder.group({
    reference: ["0", Validators.required],
    comparands: [""],
    normalize_eres: [true],
    bins: [0.2],
    bg_sub: [true],
    roi_width: [null],
    mirrored: [true],
  });

  constructor(public appData: AppData, private router: Router) { }

  close() {
    this.appData.dialogOpen.set(false);
  }

  async submit(event: Event) {
    event.preventDefault();

    const args = this.formModel.getRawValue();

    const parse = (element: string) => {
      const h = Number(element);
      if (isNaN(h)) {
        return element;
      } else {
        return h;
      }
    };

    const payload = {
      selection: parseSelection(this.appData.selection()),
      args: {
        ...args,
        reference: parse(args.reference),
        comparands: args.comparands === "" ?
          null :
          args.comparands.split(",")
            .map(element => element.trim())
            .map(element => parse(element)),
      },
    };

    const response = await fetch("http://127.0.0.1:8000/ratio_plotter", {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify(payload),
    });

    if (response.ok) {
      this.close();
    } else {
      console.error(response);
    }

    this.plot(
      await response.json(),
    );
  }

  plot(data: any) {
    const separate_detpairs = true;

    const energies = data.energies;
    const ratios = data.ratios;
    const dratios = data.dratios;

    const detpairs = Object.keys(energies);

    const n_detpairs = detpairs.length;
    const n_comparands = Object.values<Array<Array<number>>>(ratios)[0].length;
    const n_plots = separate_detpairs ? n_detpairs : 1;

    const plot_indices = separate_detpairs ?
      [...Array(n_detpairs).keys()] :
      Array(n_detpairs).fill(1);

    const plotData = Array<PlotTrace>();
    const plotAxes: PlotAxes = {};

    for (let i = 0; i < n_detpairs; i++) {
      const i_plot = plot_indices[i];
      const detpair = detpairs[i];

      plotAxes[i_plot == 0 ? "xaxis" : `xaxis${i_plot + 1}`] = {
        title: { text: "Energy - 511 keV" },
      }
      plotAxes[i_plot == 0 ? "yaxis" : `yaxis${i_plot + 1}`] = {
        title: { text: "Ratio to Reference" },
      }

      for (let j = 0; j < n_comparands; j++) {
        const trace: PlotTrace = {
          y: ratios[detpair][j] as Array<number>,
          xaxis: i_plot == 0 ? "x" : `x${i_plot + 1}`,
          yaxis: i_plot == 0 ? "y" : `y${i_plot + 1}`,
          type: "scatter",
          line: {
            shape: "spline",
            smoothing: 1.0,  // number between 0 and 1.3
          }
        };

        if (!separate_detpairs) {
          trace.name = detpair;
        }

        plotData.push(trace);
      }
    }

    this.appData.n_plots = n_plots;
    this.appData.plotAxes = plotAxes;
    this.appData.plotData.set(plotData);
    this.router.navigate(["/plots"]);
  }
}
