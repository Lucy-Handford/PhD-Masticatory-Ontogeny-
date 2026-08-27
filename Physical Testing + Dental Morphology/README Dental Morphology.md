# Dental morphology quantification workflow

This repository structure separates **application of the final method** from the
**validation and sensitivity analyses used to develop it**.

## Files

- `R/dental_morphology_workflow.R` — reusable analysis script. A new user only
  needs a landmark CSV plus corresponding upper and lower STL meshes. The final
  validated RoC/SA/CSA parameters are applied directly.
- `R/validation_and_sensitivity_analysis.R` — thesis method-development script.
  This contains repeatability analyses, Blender validation and parameter
  sensitivity testing. It is not required for routine analysis of a new specimen.
  The paths at the top are deliberately thesis-specific and should be edited to
  match whichever validation/example files are deposited with the repository.

## Suggested repository structure

```text
repository/
├── R/
│   ├── dental_morphology_workflow.R
│   └── validation_and_sensitivity_analysis.R
├── Data/
│   ├── example/                 # optional worked example for new users
│   └── validation/              # D1/A1 or equivalent thesis-validation files
└── Results/
```

## Running the reusable workflow

1. Open `R/dental_morphology_workflow.R`.
2. Edit the three paths in the **USER INPUTS** section:
   - landmark CSV
   - upper STL
   - lower STL
3. Run the script from the repository root.
4. Summary and QC CSV files are written to `Results/`.

The landmark CSV is expected to contain:
`Specimen`, `Model`, `replicate`, `side`, `tooth`, `landmark`, `x`, `y`, `z`.

## Final parameters implemented

The reusable workflow uses the parameter settings selected during the thesis
validation/sensitivity work:

- RoC cap scale: `0.15`
- RoC profile cut-off: `0.60`
- subtended-angle profile cut-off: `0.60`
- outer points per cusp side for subtended angle: `5`
- incisor CSA section scale: `0.80`
- incisor CSA bounding-box buffer: `0.10`

## Example data

The D1 dataset will be supplied as worked examples and/or as
the validation datasets accompanying the thesis.

## Important limitation

The automated incisor CSA method should be used in accordance with the
validation described in the thesis. Where adjacent teeth prevent reliable
isolation of an individual cross-section (notably in some permanent incisors),
manual extraction/validation may be required.

## Validation outputs

The validation script writes summary CSV files for tooth-dimension repeatability,
RoC/SA Blender validation, CSA Blender validation, parameter sensitivity, landmark
sensitivity, and CSA repeatability to `Results/validation/`.
