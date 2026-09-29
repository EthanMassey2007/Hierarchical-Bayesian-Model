# Suggested figure caption

**Methodological workflow for Bayesian modeling of weekly municipal dengue incidence in Rio de Janeiro, Brazil, 2017-2023.** Disease, climatic, socioeconomic, geographic, and mobility data were integrated for 92 municipalities. The temporal evaluation design included training in 2017-2021, held-out testing in 2022-2023, and rolling-origin validation. The common Bayesian fitting step, shown explicitly for all candidate models, used a negative-binomial likelihood in R-INLA. Non-spatial models (M0-M5) and spatial benchmarks (S1-S2) were extended to assess alternative spatial and mobility specifications (S3-S5), rainfall heterogeneity (S6, S8, S10, S11), and temperature heterogeneity (S7, S9). Assessment distinguished predictive performance from model fit and residual spatial structure; posterior inference characterized rainfall associations across space and time. Model groups summarize alternative specifications, rather than implying inclusion of every listed component in each model. The combined rainfall specification (S11) contains additive municipality-specific and time-varying effects, not a municipality-specific space-time interaction. Colors distinguish inputs (blue), preparation and modeling (white), and assessment and inference (teal).

Abbreviations: IDHM, municipal Human Development Index; BYM2, Besag-York-Mollie 2; INLA, integrated nested Laplace approximation; MAE, mean absolute error; RMSE, root mean squared error; WAPE, weighted absolute percentage error; DIC, deviance information criterion; WAIC, widely applicable information criterion.

# Figure specifications

- One continuous diagram, 180 x 217.8 mm, intended for full-width placement.
- Vector PDF with embedded Arial regular and bold fonts.
- Editable SVG with live text and vector shapes.
- PNG exported at 600 dpi.
- Body text approximately 8.4-9.2 pt; headings approximately 9.2 pt at 180 mm width.
- Lines approximately 0.77 pt at the supplied width.
- No meaning depends solely on color: each stage and model group is labeled.

| Function | Fill | Border |
| --- | --- | --- |
| Data inputs | Pale blue #EDF3F8 | Slate blue #56758A |
| Preparation and modeling | White #FFFFFF | Gray #737B82 |
| Assessment and inference | Pale teal #EDF5F2 | Muted teal #527B70 |

Text and arrows: charcoal #252A30. Background: white #FFFFFF.

# Scientific and submission checks

The figure summarizes the supplied manuscript; the analysis code has not been audited for this figure. Exact climate lag and geographic region taxonomy are omitted because the draft contains conflicting descriptions (six versus twelve weeks; intermediate regions versus mesoregions). Reconcile these in the manuscript before submission. The sequence summarizes analytical stages, not the implementation order of preprocessing within validation folds. Confirm how standardization, interpolation, and lag selection used training data within each fold. Verify that the named assessments match those actually performed. Final artwork dimensions and accepted file types depend on the target journal, which has not been specified. Do not reduce this full-width figure to a single column without reworking the layout.

# LaTeX insertion

Place the PDF beside the manuscript and use:

```latex
\begin{figure*}[t]
    \centering
    \includegraphics[width=\textwidth]{dengue_study_workflow.pdf}
    \caption{Insert the caption above, shortening abbreviations if defined elsewhere.}
    \label{fig:study_workflow}
\end{figure*}
```

For a one-column manuscript, use `figure` instead of `figure*`. Inspect the resulting figure size: the native width is 180 mm, and a narrower text area reduces the type size proportionally.
