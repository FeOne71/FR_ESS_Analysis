---
trigger: always_on
glob: "**/*.m"
description: Standards for RPT and Field Data Analysis
---

### Data Consistency Rules
- **Voltage Resampling**: All analysis must use linear interpolation at `0.001V` resolution.
- **Master Ruler**: Never use local peak detection for segmenting field data; always align to a global reference (`MasterRulers.mat`) for consistent feature comparison.
- **Feature Format**: Maintain the 14-column feature vector format: `[Chg_dQ(1-5), Dch_dQ(1-5), Chg_PkH, Chg_PkA, Dch_PkH, Dch_PkA]`.

### Scripting Standards
- **Sectioning**: Use `%%` for logical sections to improve readability and allow partial execution.
- **Error Handling**: Use `try-catch` blocks or conditional checks when interpolating field data near voltage window boundaries.
- **Variable Naming**: Use `pred_` prefix for estimated labels (e.g., `pred_soh`).
- **Standard Windows**:
  - Charge: [3.70 V, 3.95 V]
  - Discharge: [3.75 V, 3.88 V]
