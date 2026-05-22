import { FormControl, FormGroup, Validators } from "@angular/forms";
import { Parser } from "../../../../types";

export class ParserFormElement {
  value: string | undefined;
  key: string;
  label: string;
  required: boolean;
  order: number;
  controlType: string;
  type: string;
  options: {key: string; value: string}[];

  constructor(
    options: {
      value?: string;
      key?: string;
      label?: string;
      required?: boolean;
      order?: number;
      controlType?: string;
      type?: string;
      options?: {key: string; value: string}[];
    } = {},
  ) {
    this.value = options.value;
    this.key = options.key || '';
    this.label = options.label || '';
    this.required = !!options.required;
    this.order = options.order ?? 1;
    this.controlType = options.controlType || '';
    this.type = options.type || '';
    this.options = options.options || [];
  }
}

export class ParserNumberArg extends ParserFormElement {
  override controlType = "input";
  override type = "number";
}

export class ParserTextArg extends ParserFormElement {
  override controlType = "input";
  override type = "text";
}

export class ParserSelectArg extends ParserFormElement {
  override controlType = "select";
}

export class ParserArgContext extends ParserFormElement {
  override controlType = "span";
}

export class ParserUI {
  private elements: ParserFormElement[];
  info: string;
  name: string;
  hasRepr: boolean;
  formGroup: FormGroup;

  constructor(options: {
    name?: string,
    info?: string,
    repr?: boolean,
    elements?: ParserFormElement[],
  }) {
    this.name = options.name || "";
    this.info = options.info || "";
    this.hasRepr = !!options.repr;
    this.elements = options.elements || [];
    this.formGroup = this.toFormGroup();
  }

  display(): ParserFormElement[] {
    if (!this.name) {
      return [];
    } else if (this.hasRepr) {
      return [
        new ParserArgContext({
          value: "x ↦ ",
        }),
        ...this.elements,
      ];
    } else {
      return [
        new ParserArgContext({
          value: `x ↦ ${this.name}(x`,
        }),
        ...this.elements.flatMap(element => [
          new ParserArgContext({
            value: ", ",
          }),
          element,
        ]),
        new ParserArgContext({
          value: ")",
        }),
      ]
    }
  }

  repr(): string {
    const values = this.formGroup.getRawValue();
    if (!this.name) {
      return "";
    } else if (this.hasRepr) {
      return ["x ↦ ", ...this.elements.map(element => {
        return element.key ? values[element.key] : element.value;
      })].join("");
    } else {
      return [
        `x ↦ ${this.name}(x`,
        ...this.elements.flatMap(element => [
          ", ", element.key ? values[element.key] : element.value,
        ]),
        ")",
      ].join("");
    }
  }

  private toFormGroup(): FormGroup<any> {
    const group: any = {};

    this.elements.forEach((element) => {
      if (element.controlType == "span") return;

      group[element.key] = element.required
        ? new FormControl(element.value || '', Validators.required)
        : new FormControl(element.value || '');
    });

    return new FormGroup(group);
  }

  toPayload(): Parser {
    return {
      name: this.name,
      args: this.formGroup.getRawValue(),
      repr: this.repr(),
    };
  }
}
