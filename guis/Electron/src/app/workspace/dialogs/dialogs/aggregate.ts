import { Component, inject } from "@angular/core";
import { FormBuilder, FormArray, FormGroup, ReactiveFormsModule } from "@angular/forms";
import { Router } from "@angular/router";

import { MultiCampaign, Param, parseSelection } from "../../../types";
import { AppData } from "../../../app_data";

@Component({
  templateUrl: "aggregate.html",
  styleUrl: "../dialog.css",
  imports: [ReactiveFormsModule],
})
export abstract class AggregateDialog {
  abstract title: string;
  abstract endpoint: string;

  private FormBuilder = inject(FormBuilder);

  formModel = this.FormBuilder.group({
    params: this.FormBuilder.array([this.newParam()]),
    keep_params: this.FormBuilder.array([this.newParam()]),
  });

  private newParam() {
    return this.FormBuilder.group({
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

  get keep_params() {
    return this.formModel.get("keep_params") as FormArray;
  }

  addParam() {
    this.params.push(this.newParam());
  }

  addKeepParam() {
    this.keep_params.push(this.newParam());
  }

  constructor(public appData: AppData, protected router: Router) { }

  close() {
    this.appData.dialogOpen.set(false);
  }

  async openParserDialog(kind: number, i: number): Promise<void> {
    const parser = await this.appData.openParserDialog();
    if (parser === null) return;

    let group;
    if (kind == 0) {
      group = this.params.at(i) as FormGroup;
    } else {
      group = this.keep_params.at(i) as FormGroup;
    }

    group.patchValue({ parser: parser });
  }

  async handleResponse(response: Response) {
    const newData = await response.json() as Array<MultiCampaign>;

    this.appData.deselect.update(val => !val);
    this.appData.clearSelection();
    this.appData.data.set(newData);
  }

  async submit(event: Event): Promise<void> {
    event.preventDefault();

    const payload = {
      selection: parseSelection(this.appData.selection()),
      args: this.formModel.getRawValue(),
    };

    payload.args.keep_params = payload.args.keep_params.filter(p => p.param != "");

    const response = await fetch(this.endpoint, {
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

    this.handleResponse(response);
  }
}
