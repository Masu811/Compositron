import { Component, inject } from "@angular/core";
import { FormArray, FormBuilder, FormControl, FormGroup, NonNullableFormBuilder, ReactiveFormsModule, Validators } from "@angular/forms";
import { Param, parseSelection } from "../../../types";
import { AppData } from "../../../app_data";

@Component({
  templateUrl: "cluster.html",
  styleUrl: "../dialog.css",
  imports: [ReactiveFormsModule],
})
export class ClusterDialog {
  formBuilder = inject(FormBuilder);

  formModel = this.formBuilder.group({
    params: this.formBuilder.array([this.newParam()]),
    method: ["n_clusters"],
    clusters: this.formBuilder.array([this.formBuilder.control("")]),
    n_clusters: [null, Validators.min(0)],
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
    } as Param);
  }

  get params() {
    return this.formModel.get("params") as FormArray;
  }

  get clusters() {
    return this.formModel.get("clusters") as FormArray;
  }

  get method() {
    return this.formModel.get("method") as FormControl;
  }

  addParam() {
    this.params.push(this.newParam());
    this.clusters.push(this.formBuilder.control(""));
  }

  constructor(public appData: AppData) { }

  close() {
    this.appData.dialogOpen.set(false);
  }

  async openParserDialog(i: number): Promise<void> {
    const parser = await this.appData.openParserDialog();
    if (parser === null) return;
    const group = this.params.at(i) as FormGroup;
    group.patchValue({ parser: parser });
  }

  async submit(event: Event): Promise<void> {
    event.preventDefault();

    const args = this.formModel.getRawValue();

    const payload = {
      selection: parseSelection(this.appData.selection()),
      args: {
        ...args,
        clusters: args.clusters.map(param_cluster => {
          if (param_cluster == null || param_cluster == "") {
            return [null];
          }
          return param_cluster.split(",").map(elem => Number(elem.trim()));
        })
      }
    };

    const response = await fetch("http://127.0.0.1:8000/cluster", {
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
  }
}
