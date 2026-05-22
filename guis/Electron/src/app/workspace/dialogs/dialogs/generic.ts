import { inject } from "@angular/core";
import { FormBuilder, FormArray, FormGroup } from "@angular/forms";
import { Router } from "@angular/router";

import { Param, parseSelection } from "../../../types";
import { AppData } from "../../../app_data";

export abstract class GenericCampaignDialog {
  abstract title: string;
  abstract endpoint: string;

  private FormBuilder = inject(FormBuilder);

  formModel = this.FormBuilder.group({
    params: this.FormBuilder.array([this.newParam()])
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

  addItem() {
    this.params.push(this.newParam());
  }

  constructor(public appData: AppData, protected router: Router) { }

  close() {
    this.appData.dialogOpen.set(false);
  }

  async openParserDialog(i: number): Promise<void> {
    const parser = await this.appData.openParserDialog();
    if (parser === null) return;
    const group = this.params.at(i) as FormGroup;
    group.patchValue({ parser: parser });
  }

  abstract handleResponse(response: Response): Promise<void>;

  async submit(event: Event): Promise<void> {
    event.preventDefault();

    const payload = {
      selection: parseSelection(this.appData.selection()),
      args: this.formModel.getRawValue(),
    };

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
