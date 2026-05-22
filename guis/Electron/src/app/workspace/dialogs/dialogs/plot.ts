import { Component, inject } from "@angular/core";
import { FormArray, ReactiveFormsModule, NonNullableFormBuilder } from "@angular/forms";
import { AppData } from "../../../app_data";
import { Param, ParsedSelection, parseSelection, PlotAxes, PlotTrace } from "../../../types";
import { Router } from "@angular/router";

@Component({
  templateUrl: "plot.html",
  styleUrl: "../dialog.css",
  imports: [ReactiveFormsModule],
})
export class PlotDialog {
  private formBuilder = inject(NonNullableFormBuilder);

  formModel = this.formBuilder.group({
    x_elements: this.formBuilder.array([
      this.newParam(),
    ]),
    y_elements: this.formBuilder.array([
      this.newParam(),
    ]),
    z_elements: this.formBuilder.array([
      this.newParam(),
    ]),
    errorbars: [true],
    scatter: [false],
    filled: [false],
    label: [""],
    label_is_param: [false],
    detname_in_label: [null],
  });

  private newParam() {
    return this.formBuilder.group({
      param: "",
      target: "auto",
      parser: {
        name: "",
        args: {},
        repr: "",
      },
    });
  }

  get xParams() {
    return this.formModel.get('x_elements') as FormArray;
  }

  addxParam() {
    this.xParams.push(this.newParam());
  }

  get yParams() {
    return this.formModel.get('y_elements') as FormArray;
  }

  addyParam() {
    this.yParams.push(this.newParam());
  }

  get zParams() {
    return this.formModel.get('z_elements') as FormArray;
  }

  addzParam() {
    this.zParams.push(this.newParam());
  }

  validNumOfParams() {
    const lx = this.xParams.length;
    const ly = this.yParams.length;
    const lz = this.zParams.length;
    if (lx > 1) {
      if (ly > 1 && lx != ly) {
        return false;
      }
      if (lz > 1 && lx != lz) {
        return false;
      }
    }
    if (ly > 1) {
      if (lz > 1 && ly != lz) {
        return false;
      }
    }
    return true;
  }

  constructor(public appData: AppData, private router: Router) { }

  close() {
    this.appData.dialogOpen.set(false);
  }

  async openParserDialog(dim: string, i: number) {
    const parser = await this.appData.openParserDialog();
    if (parser === null) return;
    switch (dim) {
      case "x":
        this.xParams.at(i).patchValue({ parser: parser })
        break;
      case "y":
        this.yParams.at(i).patchValue({ parser: parser })
        break;
      case "z":
        this.zParams.at(i).patchValue({ parser: parser })
        break;
    }
  }

  async get(selection: ParsedSelection, element: Param) {
    const payload = {
      selection: selection,
      param: element,
    }

    const response = await fetch("http://127.0.0.1:8000/get", {
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
      return;
    }

    return await response.json();
  }

  async submit(event: Event) {
    event.preventDefault();

    const selection = parseSelection(this.appData.selection());
    const get = async (element: Param) => await this.get(selection, element);

    const params = this.formModel.getRawValue();

    let x_elements = params.x_elements;
    let y_elements = params.y_elements;
    let z_elements = params.z_elements;

    const nx = x_elements.length;
    const ny = y_elements.length;
    const nz = z_elements.length;

    const n_plots = Math.max(nx, ny, nz);

    if (nx < n_plots) {
      x_elements = Array(n_plots).fill(x_elements[0]);
    }
    if (ny < n_plots) {
      y_elements = Array(n_plots).fill(y_elements[0]);
    }
    if (nz < n_plots) {
      z_elements = Array(n_plots).fill(z_elements[0]);
    }

    const plotData = Array<PlotTrace>();
    const plotAxes: PlotAxes = {};

    for (let i = 0; i < n_plots; i++) {
      let y_data = await get(y_elements[i] as Param);

      const n_mc = y_data.length;

      const hasX = x_elements[i].param != "";
      const hasZ = z_elements[i].param != "";

      let x_data = [];
      let z_data = [];

      if (hasX) {
        x_data = await get(x_elements[i] as Param);
        if (x_data.length != n_mc) {
          console.error("Number of Measurement Campaigns in plot data for x and y don't match");
          return;
        }
      }
      if (hasZ) {
        z_data = await get(z_elements[i] as Param);
        if (z_data.length != n_mc) {
          console.error("Number of Measurement Campaigns in plot data for z and y don't match");
          return;
        }
      }

      plotAxes[i == 0 ? "xaxis" : `xaxis${i + 1}`] = {
        title: { text: x_elements[i]?.param || "Measurement" },
      }
      plotAxes[i == 0 ? "yaxis" : `yaxis${i + 1}`] = {
        title: { text: y_elements[i].param || "Measurement" },
      }

      for (let j = 0; j < n_mc; j++) {
        let y_labels = Object.keys(y_data[j]);
        let y_values = Object.values(y_data[j]);

        let x_labels: Array<any> = [];
        let x_values: Array<any> = [];

        let z_labels: Array<any> = [];
        let z_values: Array<any> = [];

        if (hasX) {
          x_labels = Object.keys(x_data[j]);
          x_values = Object.values(x_data[j]);
        }

        if (hasZ) {
          z_labels = Object.keys(z_data[j]);
          z_values = Object.values(z_data[j]);
        }

        const nx = x_labels.length;
        const ny = y_labels.length;
        const nz = z_labels.length;

        const n_traces = Math.max(nx, ny, nz);

        if (hasX && nx < n_traces) {
          if (ny > 1) {
            console.error("Number of traces in x data is incompatible with y and/or z");
            return;
          }
          x_labels = Array(n_traces).fill(x_labels[0]);
          x_values = Array(n_traces).fill(x_values[0]);
        }
        if (ny < n_traces) {
          if (ny > 1) {
            console.error("Number of traces in y data is incompatible with x and/or z");
            return;
          }
          y_labels = Array(n_traces).fill(y_labels[0]);
          y_values = Array(n_traces).fill(y_values[0]);
        }
        if (hasZ && nz < n_traces) {
          if (nz > 1) {
            console.error("Number of traces in z data is incompatible with x and/or y");
            return;
          }
          z_labels = Array(n_traces).fill(z_labels[0]);
          z_values = Array(n_traces).fill(z_values[0]);
        }

        for (let k = 0; k < n_traces; k++) {
          const trace: PlotTrace = {
            y: y_values[k] as Array<number>,
            xaxis: i == 0 ? "x" : `x${i + 1}`,
            yaxis: i == 0 ? "y" : `y${i + 1}`,
          };

          if (hasX) {
            trace.x = x_values[k] as Array<number>;
          }

          if (hasZ) {
            trace.marker = {
              color: z_values[k] as Array<number>,
            };
          }

          if (ny > 1 && nx <= 1 && nz <= 1) {
            trace.name = y_labels[k];
          } else if (ny == 1 && nx > 1 && nz <= 1) {
            trace.name = x_labels[k];
          } else if (ny == 1 && nx <= 1 && nz > 1) {
            trace.name = z_labels[k];
          } else if (ny > 1 && nx > 1 && nz <= 1) {
            trace.name = `${x_labels[k]} x ${y_labels[k]}`;
          } else if (ny > 1 && nx <= 1 && nz > 1) {
            trace.name = `${y_labels[k]} x ${z_labels[k]}`;
          } else if (ny == 1 && nx > 1 && nz > 1) {
            trace.name = `${x_labels[k]} x ${z_labels[k]}`;
          }

          plotData.push(trace);
        }
      }
    }

    this.appData.n_plots = n_plots;
    this.appData.plotAxes = plotAxes;
    this.appData.plotData.set(plotData);
    this.router.navigate(["/plots"]);
  }
}
