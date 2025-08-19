%%{init: {'theme':'default', 'themeVariables': { 'primaryColor':'#fff', 'primaryTextColor':'#000', 'primaryBorderColor':'#7C0000', 'lineColor':'#000', 'secondaryColor':'#f9f9f9', 'tertiaryColor':'#fff', 'background':'#ffffff', 'mainBkg':'#e8e8e8', 'secondBkg':'#f0f0f0', 'tertiaryBkg':'#ffffff'}}}%%
graph TD
    A[Raw RNA Reads] --> B[Quality Filtering]
    B --> C{Analysis Type}
    C -->|AIV Analysis| D[Map to AIV References]
    C -->|Virome Analysis| E[DIAMOND BLASTx]
    D --> F[Select Best Reference]
    F --> G[Variant Calling]
    G --> H[Consensus Generation]
    E --> I[Viral Family Classification]
    I --> J[Abundance Heatmap]
    H --> K[Complete AIV Genome]
    J --> L[Virome Profile]
