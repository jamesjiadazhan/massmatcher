`%||%` <- function(x, y) {
    if (is.null(x) || length(x) == 0) {
        return(y)
    }
    x
}

locate_massmatcher_app_dir <- function() {
    current_file <- tryCatch(
        normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
        error = function(e) ""
    )

    app_dir_candidates <- unique(c(
        if (nzchar(current_file)) dirname(current_file),
        normalizePath(file.path(getwd(), "inst", "shiny", "massmatcher-app"), winslash = "/", mustWork = FALSE),
        normalizePath(file.path(getwd(), "massmatcher", "inst", "shiny", "massmatcher-app"), winslash = "/", mustWork = FALSE)
    ))
    app_dir_candidates <- app_dir_candidates[dir.exists(app_dir_candidates)]
    app_dir_candidates[1]
}

locate_massmatcher_source_root <- function() {
    current_file <- tryCatch(
        normalizePath(sys.frame(1)$ofile, winslash = "/", mustWork = TRUE),
        error = function(e) ""
    )

    app_dir_candidates <- unique(c(
        if (nzchar(current_file)) dirname(current_file),
        normalizePath(file.path(getwd(), "inst", "shiny", "massmatcher-app"), winslash = "/", mustWork = FALSE),
        normalizePath(file.path(getwd(), "massmatcher", "inst", "shiny", "massmatcher-app"), winslash = "/", mustWork = FALSE),
        normalizePath(getwd(), winslash = "/", mustWork = FALSE)
    ))
    app_dir_candidates <- app_dir_candidates[dir.exists(app_dir_candidates)]

    root_candidates <- unique(unlist(lapply(app_dir_candidates, function(app_dir) {
        c(
            normalizePath(file.path(app_dir, "..", "..", ".."), winslash = "/", mustWork = FALSE),
            normalizePath(file.path(app_dir, "..", "..", "..", "massmatcher"), winslash = "/", mustWork = FALSE),
            normalizePath(app_dir, winslash = "/", mustWork = FALSE)
        )
    }), use.names = FALSE))
    root_candidates <- root_candidates[dir.exists(root_candidates)]

    for (root in root_candidates) {
        r_dir <- file.path(root, "R")
        if (file.exists(file.path(root, "DESCRIPTION")) &&
            dir.exists(file.path(root, "inst")) &&
            dir.exists(file.path(root, "data")) &&
            dir.exists(r_dir) &&
            length(list.files(r_dir, pattern = "\\.[Rr]$", full.names = TRUE)) > 0) {
            return(root)
        }
    }

    NULL
}

bootstrap_massmatcher_backend <- function() {
    source_root <- locate_massmatcher_source_root()

    if (!is.null(source_root)) {
        options(massmatcher.bundle_root = source_root)

        import_env <- new.env(parent = baseenv())
        import_packages <- c(
            "data.table",
            "MetaboCoreUtilsAdduct",
            "arrow",
            "classyfireR",
            "purrr",
            "readr",
            "preprocessCore",
            "imputeLCMD",
            "magrittr",
            "rlang",
            "tibble",
            "tidyr",
            "dplyr"
        )
        for (pkg in import_packages) {
            if (!requireNamespace(pkg, quietly = TRUE)) {
                next
            }
            pkg_exports <- getNamespaceExports(pkg)
            for (export_name in pkg_exports) {
                assign(export_name, getExportedValue(pkg, export_name), envir = import_env)
            }
        }

        pkg_env <- new.env(parent = import_env)
        r_files <- sort(list.files(file.path(source_root, "R"), pattern = "\\.[Rr]$", full.names = TRUE))
        r_files <- r_files[!grepl("/run_massmatcher_app\\.R$", r_files)]
        for (r_file in r_files) {
            sys.source(r_file, envir = pkg_env)
        }

        return(list(mode = "bundle", root = source_root, env = pkg_env))
    }

    if (requireNamespace("massmatcher", quietly = TRUE)) {
        return(list(mode = "package", root = NULL, env = NULL))
    }

    stop("Could not locate the massmatcher source bundle or an installed massmatcher package.")
}

massmatcher_backend <- bootstrap_massmatcher_backend()

mm_call <- function(name, ...) {
    if (identical(massmatcher_backend$mode, "bundle")) {
        return(get(name, envir = massmatcher_backend$env, inherits = FALSE)(...))
    }

    getExportedValue("massmatcher", name)(...)
}

DEMO_MZ_MASS_MATCH <- c(160.1332, 176.1282)
DEMO_MZ_ISOTOPOLOGUE <- c(161.1366, 177.1315)
DEMO_TARGET_MASS <- c(180.063388, 202.045322)

parse_adducts <- function(x, default = "M+H", allow_null = FALSE) {
    if (is.null(x) || length(x) == 0) {
        if (isTRUE(allow_null)) {
            return(NULL)
        }
        return(default)
    }

    if (length(x) > 1) {
        values <- trimws(as.character(x))
        values <- values[values != ""]
        if (length(values) == 0) {
            if (isTRUE(allow_null)) {
                return(NULL)
            }
            return(default)
        }
        return(unique(values))
    }

    x <- trimws(as.character(x))
    if (identical(x, "")) {
        if (isTRUE(allow_null)) {
            return(NULL)
        }
        return(default)
    }

    values <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
    values <- trimws(values)
    values <- values[values != ""]
    unique(values)
}

default_search_adduct <- function(ion_mode = c("positive", "negative")) {
    ion_mode <- match.arg(ion_mode)
    if (identical(ion_mode, "negative")) "M-H" else "M+H"
}

default_cluster_adducts <- function(ion_mode = c("positive", "negative")) {
    ion_mode <- match.arg(ion_mode)
    if (identical(ion_mode, "negative")) {
        return(c("M-H", "M+Cl", "M-H-H2O", "2M-H"))
    }

    c(
        "M+H", "M+Na", "M+2Na-H", "M+H-H2O", "M+H-NH3",
        "M+ACN+H", "M+ACN+2H", "2M+H", "M+2H", "M+H-2H2O"
    )
}

get_supported_search_adducts <- function(ion_mode = c("positive", "negative")) {
    ion_mode <- match.arg(ion_mode)
    adduct_env <- new.env(parent = emptyenv())
    utils::data("adduct_definition", package = "MetaboCoreUtilsAdduct", envir = adduct_env)
    adduct_definition <- get("adduct_definition", envir = adduct_env, inherits = FALSE)

    preferred <- if (identical(ion_mode, "positive")) {
        c(
            "M+H", "M+Na", "M+K", "M+NH4", "M+2H",
            "M+ACN+H", "M+ACN+2H", "M+H-H2O", "M+H-NH3",
            "2M+H", "M+2Na-H", "M+2K-H"
        )
    } else {
        c("M-H", "M+Cl", "M+FA-H", "M+Hac-H", "2M-H", "M-H-H2O")
    }

    available <- intersect(preferred, adduct_definition$name)
    if (length(available) == 0) {
        available <- unique(adduct_definition$name)
    }
    available
}

parse_numeric_values <- function(x) {
    x <- trimws(x %||% "")
    if (identical(x, "")) {
        stop("Please provide one or more numeric values.")
    }

    values <- unlist(strsplit(x, "[[:space:],;]+"), use.names = FALSE)
    values <- trimws(values)
    values <- values[values != ""]
    numeric_values <- suppressWarnings(as.numeric(values))

    if (length(numeric_values) == 0 || any(is.na(numeric_values))) {
        stop("Values must be numeric (one per line is recommended).")
    }

    numeric_values
}

extract_mz_values_from_table <- function(table, context_label = "uploaded table") {
    table <- as.data.frame(table)
    if (ncol(table) < 1) {
        stop(context_label, " must include at least one column with m/z values.")
    }

    mz_values <- suppressWarnings(as.numeric(table[[1]]))
    mz_values <- mz_values[!is.na(mz_values)]
    if (length(mz_values) == 0) {
        stop("No numeric m/z values found in the first column of ", context_label, ".")
    }

    mz_values
}

read_mz_values_from_upload <- function(uploaded_file, context_label, require_single_column = TRUE) {
    uploaded <- read_uploaded_table(uploaded_file)
    if (isTRUE(require_single_column) && ncol(as.data.frame(uploaded)) != 1) {
        stop(context_label, " must have exactly one column named `mz`.")
    }

    extract_mz_values_from_table(table = uploaded, context_label = context_label)
}

mz_values_to_unknown_feature <- function(mz_values) {
    data.frame(
        mz = as.numeric(mz_values),
        time = seq_along(mz_values),
        stringsAsFactors = FALSE
    )
}

mz_values_to_lines <- function(mz_values) {
    paste(format(as.numeric(mz_values), scientific = FALSE, trim = TRUE), collapse = "\n")
}

read_uploaded_table <- function(uploaded_file) {
    if (is.null(uploaded_file)) {
        stop("Please upload a file.")
    }

    extension <- tolower(tools::file_ext(uploaded_file$name))

    if (extension == "csv") {
        return(readr::read_csv(uploaded_file$datapath, show_col_types = FALSE, progress = FALSE))
    }
    if (extension %in% c("tsv", "txt")) {
        return(readr::read_tsv(uploaded_file$datapath, show_col_types = FALSE, progress = FALSE))
    }
    if (extension == "rds") {
        return(readRDS(uploaded_file$datapath))
    }

    stop("Unsupported file format. Please upload CSV, TSV/TXT, or RDS.")
}

coerce_mz_time_columns <- function(table) {
    table <- as.data.frame(table)
    if (ncol(table) < 2) {
        stop("Input table must have at least two columns (mz and time).")
    }

    colnames(table)[1:2] <- c("mz", "time")
    table$mz <- suppressWarnings(as.numeric(table$mz))
    table$time <- suppressWarnings(as.numeric(table$time))

    if (all(is.na(table$mz))) {
        stop("The first column (mz) could not be parsed as numeric values.")
    }
    if (all(is.na(table$time))) {
        stop("The second column (time) could not be parsed as numeric values.")
    }

    table <- table[!is.na(table$mz) & !is.na(table$time), , drop = FALSE]
    if (nrow(table) == 0) {
        stop("No valid rows remained after parsing mz/time as numeric.")
    }

    table
}

load_demo_feature_table <- function(rows = 500, cols = 12) {
    env <- new.env(parent = emptyenv())

    if (identical(massmatcher_backend$mode, "bundle")) {
        data_file <- file.path(massmatcher_backend$root, "data", "feature_table_exp_hilicpos.rda")
        if (!file.exists(data_file)) {
            stop("Bundled demo data `feature_table_exp_hilicpos` was not found in the deployed source bundle.")
        }
        load(data_file, envir = env)
    } else {
        data("feature_table_exp_hilicpos", package = "massmatcher", envir = env)
    }

    if (!exists("feature_table_exp_hilicpos", envir = env, inherits = FALSE)) {
        stop("Bundled demo data `feature_table_exp_hilicpos` was not found.")
    }

    demo <- get("feature_table_exp_hilicpos", envir = env, inherits = FALSE)
    demo <- as.data.frame(demo)
    demo <- demo[seq_len(min(rows, nrow(demo))), seq_len(min(cols, ncol(demo))), drop = FALSE]
    colnames(demo)[1:2] <- c("mz", "time")
    demo
}

make_datatable <- function(table) {
    DT::datatable(
        table,
        rownames = FALSE,
        filter = "top",
        options = list(
            pageLength = 20,
            scrollX = TRUE,
            lengthMenu = c(20, 50, 100)
        )
    )
}

make_grouped_search_datatable <- function(table, page_size = 10) {
    page_size <- as.integer(page_size)
    if (!is.finite(page_size) || !(page_size %in% c(10L, 25L, 50L, 100L))) {
        page_size <- 10L
    }

    widget <- DT::datatable(
        table,
        filter = "none",
        rownames = FALSE,
        fillContainer = FALSE,
        width = "100%",
        options = list(
            pageLength = page_size,
            lengthMenu = list(c(10, 25, 50, 100), c("10", "25", "50", "100")),
            paging = TRUE,
            searching = TRUE,
            info = TRUE,
            scrollX = TRUE,
            dom = "lfrtip",
            ordering = TRUE,
            autoWidth = TRUE
        )
    )

    htmltools::bindFillRole(widget, item = FALSE, container = FALSE, overwrite = TRUE)
}

prepare_export_table <- function(table) {
    df <- as.data.frame(table, stringsAsFactors = FALSE)
    list_cols <- vapply(df, is.list, logical(1))
    if (any(list_cols)) {
        df[list_cols] <- lapply(df[list_cols], function(col) {
            vapply(col, function(x) {
                if (is.null(x) || length(x) == 0) {
                    return(NA_character_)
                }
                paste(as.character(unlist(x, use.names = FALSE)), collapse = ";")
            }, character(1))
        })
    }
    df
}

ui <- bslib::page_navbar(
    title = shiny::tags$div(
        style = "display:flex;align-items:center;gap:10px;",
        shiny::tags$span("massmatcher")
    ),
    theme = bslib::bs_theme(version = 5, bootswatch = "flatly", primary = "#50C878"),
    header = shiny::tagList(
        shiny::tags$style(shiny::HTML("
          .nav-pills .nav-link.active,
          .navbar .nav-link.active {
            color: inherit !important;
            background-color: transparent !important;
          }
          .search-group-card {
            border: 1px solid #d9d9d9;
            border-radius: 12px;
            background: #ffffff;
            margin-bottom: 1rem;
            overflow: hidden;
          }
          .search-group-header {
            padding: 14px 18px;
            border-bottom: 1px solid #e7e7e7;
            background: #fbfbfb;
            font-size: 1.1rem;
          }
          .search-group-body {
            padding: 14px 18px 6px 18px;
            display: block !important;
            overflow: visible;
          }
          .search-group-body > .html-widget,
          .search-group-body > .html-fill-item,
          .search-group-body > .datatables {
            display: block !important;
            flex: 0 0 auto !important;
            height: auto !important;
            min-height: 0 !important;
            max-height: none !important;
          }
          .search-group-body .dataTables_wrapper {
            width: 100%;
            height: auto !important;
          }
          .search-group-body .datatables,
          .search-group-body .html-widget,
          .search-group-body table.dataTable {
            width: 100% !important;
          }
          .search-group-body .dataTables_scrollBody {
            overflow: auto !important;
          }
          .metadata-progress-panel {
            display: none;
            margin-bottom: 1rem;
            padding: 0.85rem 1rem 0.95rem 1rem;
            border: 1px solid #d9ece0;
            border-radius: 12px;
            background: #f4fcf7;
          }
          .metadata-progress-label {
            font-weight: 600;
            margin-bottom: 0.35rem;
            color: #2f5d46;
          }
          .metadata-progress-detail {
            font-size: 0.95rem;
            color: #4c6a5a;
            margin-bottom: 0.6rem;
          }
          .metadata-progress-panel .progress {
            height: 0.95rem;
            background-color: #dff4e7;
          }
          .metadata-progress-panel .progress-bar {
            background-color: #50C878;
          }
        ")),
        shiny::tags$script(shiny::HTML("
          Shiny.addCustomMessageHandler('massmatcher-metadata-progress', function(message) {
            var panel = document.getElementById('metadata-progress-panel');
            var label = document.getElementById('metadata-progress-label');
            var detail = document.getElementById('metadata-progress-detail');
            var bar = document.getElementById('metadata-progress-bar');
            if (!panel || !label || !detail || !bar) return;

            var value = Number(message.value);
            if (!Number.isFinite(value)) value = 0;
            value = Math.max(0, Math.min(100, value));

            panel.style.display = message.visible === false ? 'none' : 'block';
            label.textContent = message.label || 'Additional Meta Data progress';
            detail.textContent = message.detail || '';
            bar.style.width = value + '%';
            bar.setAttribute('aria-valuenow', value);
            bar.textContent = value + '%';

            bar.classList.toggle('progress-bar-striped', !!message.striped);
            bar.classList.toggle('progress-bar-animated', !!message.animated);
          });
        "))
    ),
    bslib::nav_panel(
        title = "Search",
        bslib::layout_sidebar(
            sidebar = bslib::sidebar(
                id = "search_sidebar",
                open = TRUE,
                width = 360,
                shiny::radioButtons(
                    inputId = "search_mode",
                    label = "Search method",
                    choices = c(
                        "m/z matching" = "mass_match",
                        "Isotopologue m/z matching" = "isotopologue_mass_match",
                        "Target monoisotopic mass search" = "target_mass_search"
                    ),
                    selected = "mass_match"
                ),
                shiny::selectInput(
                    inputId = "search_database",
                    label = "Database",
                    choices = c("metorigindb", "pubchem"),
                    selected = "metorigindb"
                ),
                shiny::selectInput(
                    inputId = "search_ion_mode",
                    label = "Ion mode",
                    choices = c("positive", "negative"),
                    selected = "positive"
                ),
                shiny::selectizeInput(
                    inputId = "search_adducts",
                    label = "Adduct(s)",
                    choices = get_supported_search_adducts("positive"),
                    selected = default_search_adduct("positive"),
                    multiple = TRUE,
                    options = list(plugins = list("remove_button"))
                ),
                shiny::numericInput(
                    inputId = "search_ppm",
                    label = "m/z ppm threshold",
                    value = 10,
                    min = 0.1,
                    step = 0.1
                ),
                shiny::conditionalPanel(
                    condition = "input.search_mode === 'target_mass_search'",
                    shiny::textAreaInput(
                        inputId = "search_target_masses",
                        label = "Target mass(es), one per line",
                        value = "180.063388\n202.045322",
                        rows = 6
                    ),
                    shiny::fileInput(
                        inputId = "search_target_feature_file",
                        label = "Feature table upload required (CSV/TSV/RDS; first column mz, second column time)"
                    ),
                    shiny::helpText("Target monoisotopic mass search matches monoisotopic mass against m/z values within your uploaded feature table.")
                ),
                shiny::conditionalPanel(
                    condition = "input.search_mode !== 'target_mass_search'",
                    shiny::textAreaInput(
                        inputId = "search_mz_values",
                        label = "Unknown feature m/z values, one per line",
                        value = "",
                        rows = 10
                    ),
                    shiny::fileInput(
                        inputId = "search_unknown_file",
                        label = "Optional upload (CSV/TSV/RDS) with only 1 column: mz"
                    ),
                    shiny::helpText("For upload, use one column only (m/z).")
                ),
                shiny::actionButton("search_use_demo", "Load demo m/z into text box"),
                shiny::actionButton("run_search", "Run Search", class = "btn-primary")
            ),
            bslib::card(
                bslib::card_header(
                    shiny::tags$div(
                        style = "display:flex;justify-content:space-between;align-items:center;width:100%;gap:10px;",
                        shiny::tags$div(
                            style = "display:flex;align-items:center;gap:10px;",
                            shiny::actionButton(
                                inputId = "toggle_search_sidebar",
                                label = "Search options",
                                icon = shiny::icon("sliders"),
                                class = "btn-outline-secondary btn-sm"
                            ),
                            shiny::tags$span("Search results")
                        ),
                        shiny::uiOutput("search_download_ui")
                    )
                ),
                bslib::card_body(
                    fill = FALSE,
                    shiny::textOutput("search_status"),
                    shiny::uiOutput("search_group_page_size_ui"),
                    shiny::uiOutput("search_grouped_ui")
                )
            )
        )
    ),
    bslib::nav_panel(
        title = "Clustering",
        bslib::layout_sidebar(
            sidebar = bslib::sidebar(
                width = 380,
                shiny::fileInput(
                    inputId = "cluster_feature_file",
                    label = "Feature table for clustering (CSV/TSV/RDS)"
                ),
                shiny::actionButton("cluster_use_demo", "Use bundled demo features"),
                shiny::uiOutput("cluster_feature_source_ui"),
                shiny::selectInput(
                    inputId = "cluster_database",
                    label = "Database",
                    choices = c("metorigindb", "pubchem"),
                    selected = "metorigindb"
                ),
                shiny::numericInput(
                    inputId = "cluster_mz_threshold",
                    label = "m/z threshold (ppm)",
                    value = 10,
                    min = 0.1,
                    step = 0.1
                ),
                shiny::selectInput(
                    inputId = "cluster_ion_mode",
                    label = "Ion mode",
                    choices = c("positive", "negative"),
                    selected = "positive"
                ),
                shiny::selectizeInput(
                    inputId = "cluster_adducts",
                    label = "Adduct(s)",
                    choices = get_supported_search_adducts("positive"),
                    selected = default_cluster_adducts("positive"),
                    multiple = TRUE,
                    options = list(plugins = list("remove_button"))
                ),
                shiny::selectInput(
                    inputId = "cluster_imputation",
                    label = "Imputation",
                    choices = c("half_min", "QRILC", "None"),
                    selected = "half_min"
                ),
                shiny::numericInput(
                    inputId = "cluster_adduct_r",
                    label = "Adduct correlation r threshold",
                    value = 0.39,
                    min = -1,
                    max = 1,
                    step = 0.01
                ),
                shiny::numericInput(
                    inputId = "cluster_adduct_time",
                    label = "Adduct RT window (s)",
                    value = 6,
                    min = 0,
                    step = 0.5
                ),
                shiny::numericInput(
                    inputId = "cluster_isotope_r",
                    label = "Isotope correlation r threshold",
                    value = 0.71,
                    min = -1,
                    max = 1,
                    step = 0.01
                ),
                shiny::numericInput(
                    inputId = "cluster_isotope_time",
                    label = "Isotope RT window (s)",
                    value = 4,
                    min = 0,
                    step = 0.5
                ),
                shiny::actionButton("run_cluster", "Run Clustering", class = "btn-primary"),
                shiny::hr(),
                shiny::selectInput(
                    inputId = "cluster_view_table",
                    label = "Result table",
                    choices = c("mz_only", "mz_only_isotope", "clustering_all", "clustering_main"),
                    selected = "clustering_main"
                )
            ),
            bslib::card(
                bslib::card_header(
                    shiny::tags$div(
                        style = "display:flex;justify-content:space-between;align-items:center;width:100%;gap:10px;",
                        shiny::tags$span("Clustering results"),
                        shiny::uiOutput("cluster_download_ui")
                    )
                ),
                bslib::card_body(
                    shiny::textOutput("cluster_status"),
                    DT::DTOutput("cluster_table")
                )
            )
        )
    ),
    bslib::nav_panel(
        title = "Additional Meta Data",
        bslib::layout_sidebar(
            sidebar = bslib::sidebar(
                width = 360,
                shiny::fileInput(
                    inputId = "meta_upload_file",
                    label = "Upload data to add meta data (CSV/TSV/RDS)"
                ),
                shiny::helpText("This module adds HMDB metabolite concentration for the selected biospecimen and ClassyFire chemical classification information."),
                shiny::helpText("HMDB concentration values represent the average of all reported normal metabolite concentrations across all populations in the packaged HMDB dataset."),
                shiny::helpText("Upload a table containing InChIKey and/or HMDB_ID columns."),
                shiny::selectInput(
                    inputId = "meta_database",
                    label = "Database",
                    choices = c("metorigindb", "pubchem"),
                    selected = "metorigindb"
                ),
                shiny::selectInput(
                    inputId = "meta_biospecimen",
                    label = "Biospecimen",
                    choices = c("Blood"),
                    selected = "Blood"
                ),
                shiny::actionButton("run_metadata", "Add Meta Data", class = "btn-primary")
            ),
            bslib::card(
                bslib::card_header(
                    shiny::tags$div(
                        style = "display:flex;justify-content:space-between;align-items:center;width:100%;gap:10px;",
                        shiny::tags$span("Additional metadata results"),
                        shiny::uiOutput("metadata_download_ui")
                    )
                ),
                bslib::card_body(
                    shiny::tags$div(
                        id = "metadata-progress-panel",
                        class = "metadata-progress-panel",
                        shiny::tags$div(
                            id = "metadata-progress-label",
                            class = "metadata-progress-label",
                            "Additional Meta Data progress"
                        ),
                        shiny::tags$div(
                            id = "metadata-progress-detail",
                            class = "metadata-progress-detail",
                            "Waiting to start."
                        ),
                        shiny::tags$div(
                            class = "progress",
                            shiny::tags$div(
                                id = "metadata-progress-bar",
                                class = "progress-bar progress-bar-striped progress-bar-animated",
                                role = "progressbar",
                                style = "width:0%;",
                                `aria-valuemin` = "0",
                                `aria-valuemax` = "100",
                                `aria-valuenow` = "0",
                                "0%"
                            )
                        )
                    ),
                    shiny::textOutput("metadata_status"),
                    DT::DTOutput("metadata_table")
                )
            )
        )
    )
)

server <- function(input, output, session) {
    state <- shiny::reactiveValues(
        search_result = NULL,
        search_mode_used = NULL,
        search_grouping = NULL,
        search_run_stamp = NULL,
        clustering_result = NULL,
        metadata_table = NULL,
        cluster_use_demo = FALSE,
        cluster_demo_table = NULL,
        search_status = "Ready.",
        cluster_status = "Ready.",
        metadata_status = "Ready. Add biospecimen-specific metabolite concentration and ClassyFire chemical classification to an uploaded table."
    )

    output$search_status <- shiny::renderText(state$search_status)
    output$cluster_status <- shiny::renderText(state$cluster_status)
    output$metadata_status <- shiny::renderText(state$metadata_status)

    update_metadata_progress <- function(value, detail, label = "Additional Meta Data progress", visible = TRUE, striped = TRUE, animated = TRUE) {
        numeric_value <- suppressWarnings(as.numeric(value))
        if (!is.finite(numeric_value)) {
            numeric_value <- 0
        }
        numeric_value <- max(0, min(100, numeric_value))

        session$sendCustomMessage(
            "massmatcher-metadata-progress",
            list(
                value = numeric_value,
                detail = detail %||% "",
                label = label,
                visible = isTRUE(visible),
                striped = isTRUE(striped),
                animated = isTRUE(animated)
            )
        )
    }

    output$search_download_ui <- shiny::renderUI({
        shiny::req(state$search_result)
        shiny::downloadButton("download_search", "Download search table")
    })

    output$cluster_download_ui <- shiny::renderUI({
        shiny::req(state$clustering_result)
        shiny::tagList(
            shiny::downloadButton("download_cluster", "Download selected table"),
            shiny::downloadButton("download_cluster_all", "Download all (ZIP)")
        )
    })

    output$metadata_download_ui <- shiny::renderUI({
        shiny::req(state$metadata_table)
        shiny::downloadButton("download_metadata", "Download metadata table")
    })

    output$cluster_feature_source_ui <- shiny::renderUI({
        if (isTRUE(state$cluster_use_demo)) {
            return(
                shiny::tagList(
                    shiny::tags$small("Using bundled demo feature table: feature_table_exp_hilicpos (500 x 12 subset)."),
                    shiny::tags$br(),
                    shiny::downloadLink("download_cluster_demo", "Download demo feature table (CSV)")
                )
            )
        }

        if (!is.null(input$cluster_feature_file)) {
            return(
                shiny::tags$small(paste0("Using uploaded feature table: ", input$cluster_feature_file$name))
            )
        }

        shiny::tags$small("No clustering feature table selected yet.")
    })

    available_metadata_biospecimens <- function() {
        biospecimen_choices <- tryCatch(
            {
                concentration_df <- mm_call("hmdb_concentration_loading")
                if (!"Biospecimen" %in% colnames(concentration_df)) {
                    return("Blood")
                }
                choices <- sort(unique(as.character(concentration_df$Biospecimen)))
                choices <- choices[!is.na(choices) & choices != ""]
                if (length(choices) == 0) {
                    return("Blood")
                }
                if ("Blood" %in% choices) {
                    choices <- c("Blood", setdiff(choices, "Blood"))
                }
                choices
            },
            error = function(e) {
                "Blood"
            }
        )

        selected_choice <- if ("Blood" %in% biospecimen_choices) "Blood" else biospecimen_choices[[1]]
        list(choices = biospecimen_choices, selected = selected_choice)
    }

    shiny::observe({
        biospecimen_config <- available_metadata_biospecimens()
        shiny::updateSelectInput(
            session = session,
            inputId = "meta_biospecimen",
            choices = biospecimen_config$choices,
            selected = biospecimen_config$selected
        )
    })

    shiny::observeEvent(input$search_ion_mode, {
        mode <- input$search_ion_mode %||% "positive"
        adduct_choices <- get_supported_search_adducts(mode)
        default_choice <- default_search_adduct(mode)
        selected <- if (default_choice %in% adduct_choices) default_choice else adduct_choices[[1]]
        shiny::updateSelectizeInput(
            session = session,
            inputId = "search_adducts",
            choices = adduct_choices,
            selected = selected,
            server = TRUE
        )
    }, ignoreInit = TRUE)

    shiny::observeEvent(input$cluster_ion_mode, {
        mode <- input$cluster_ion_mode %||% "positive"
        adduct_choices <- get_supported_search_adducts(mode)
        default_choices <- intersect(default_cluster_adducts(mode), adduct_choices)
        if (length(default_choices) == 0) {
            default_choices <- adduct_choices
        }

        shiny::updateSelectizeInput(
            session = session,
            inputId = "cluster_adducts",
            choices = adduct_choices,
            selected = default_choices,
            server = TRUE
        )
    }, ignoreInit = TRUE)

    parse_text_or_upload_mz <- function(text_value, uploaded_file, context_label) {
        if (!is.null(uploaded_file)) {
            return(read_mz_values_from_upload(uploaded_file, context_label = context_label))
        }

        text_value <- trimws(text_value %||% "")
        if (!identical(text_value, "")) {
            return(parse_numeric_values(text_value))
        }

        stop(
            "Please enter m/z values (one per line) or upload a file with one column: mz."
        )
    }

    should_use_search_broad_prefilter <- function(mz_values, adducts, uploaded_file = NULL) {
        if (!is.null(uploaded_file)) {
            return(TRUE)
        }

        mz_count <- length(mz_values %||% numeric(0))
        adduct_count <- length(adducts %||% character(0))

        isTRUE(mz_count > 10L) || isTRUE(adduct_count > 1L)
    }

    detect_search_group_column <- function(df, mode = NULL) {
        if (identical(mode, "target_mass_search") && "mass" %in% names(df)) {
            return("mass")
        }
        if ("mz" %in% names(df)) {
            return("mz")
        }
        if ("mass" %in% names(df)) {
            return("mass")
        }
        if ("theoretical_mz" %in% names(df)) {
            return("theoretical_mz")
        }
        NULL
    }

    format_group_value <- function(x, digits = 6) {
        x_num <- suppressWarnings(as.numeric(x))
        if (is.finite(x_num)) {
            return(format(round(x_num, digits), scientific = FALSE, trim = TRUE))
        }
        as.character(x)
    }

    search_group_page_size <- shiny::reactive({
        page_size <- suppressWarnings(as.integer(input$search_group_page_size %||% "10"))
        if (!is.finite(page_size) || !(page_size %in% c(10L, 25L, 50L, 100L))) {
            return(10L)
        }
        page_size
    })

    should_group_search_results <- function(mode_used = NULL) {
        mode_used <- mode_used %||% state$search_mode_used %||% input$search_mode
        grouping_info <- state$search_grouping

        if (!mode_used %in% c("mass_match", "isotopologue_mass_match")) {
            return(TRUE)
        }

        if (is.null(grouping_info)) {
            return(TRUE)
        }

        if (isTRUE(grouping_info$uploaded)) {
            return(FALSE)
        }

        query_count <- suppressWarnings(as.integer(grouping_info$query_count %||% NA_integer_))
        if (is.finite(query_count) && query_count > 10L) {
            return(FALSE)
        }

        TRUE
    }

    build_search_group_specs <- function() {
        shiny::req(state$search_result)
        df <- as.data.frame(state$search_result, stringsAsFactors = FALSE)

        if (nrow(df) == 0) {
            return(list(type = "empty"))
        }

        mode_used <- state$search_mode_used %||% input$search_mode
        group_col <- detect_search_group_column(df, mode = mode_used)

        if (is.null(group_col) || !should_group_search_results(mode_used)) {
            return(list(type = "single"))
        }

        group_values <- unique(df[[group_col]])
        group_values_num <- suppressWarnings(as.numeric(group_values))
        if (length(group_values) > 1 && all(is.finite(group_values_num))) {
            group_values <- group_values[order(group_values_num)]
        }

        groups <- lapply(seq_along(group_values), function(i) {
            group_value <- group_values[[i]]
            section_df <- df[df[[group_col]] == group_value, , drop = FALSE]

            if ("mz_ppm" %in% names(section_df)) {
                section_df <- section_df[order(section_df$mz_ppm), , drop = FALSE]
            }

            heading <- if (identical(group_col, "mass")) {
                paste0("MS search for target monoisotopic mass ", format_group_value(group_value))
            } else {
                paste0("MS search for ", format_group_value(group_value), " m/z")
            }

            list(
                heading = heading,
                data = prepare_export_table(section_df)
            )
        })

        list(type = "grouped", groups = groups)
    }

    search_group_specs <- shiny::reactive({
        build_search_group_specs()
    })

    shiny::observeEvent(input$search_use_demo, {
        if (identical(input$search_mode, "target_mass_search")) {
            shiny::updateTextAreaInput(session, "search_target_masses", value = mz_values_to_lines(DEMO_TARGET_MASS))
            state$search_status <- "Loaded demo target masses. Please upload a feature table (first column mz, second column time)."
        } else if (identical(input$search_mode, "isotopologue_mass_match")) {
            shiny::updateTextAreaInput(session, "search_mz_values", value = mz_values_to_lines(DEMO_MZ_ISOTOPOLOGUE))
            state$search_status <- paste0("Loaded ", length(DEMO_MZ_ISOTOPOLOGUE), " demo m/z values (same as README isotopologue example).")
        } else {
            shiny::updateTextAreaInput(session, "search_mz_values", value = mz_values_to_lines(DEMO_MZ_MASS_MATCH))
            state$search_status <- paste0("Loaded ", length(DEMO_MZ_MASS_MATCH), " demo m/z values (same as README mass_match example).")
        }
    })

    shiny::observeEvent(input$search_unknown_file, {
        req(input$search_unknown_file)
        try({
            mz_values <- read_mz_values_from_upload(
                input$search_unknown_file,
                context_label = "search upload",
                require_single_column = TRUE
            )
            state$search_status <- paste0(
                "Uploaded file detected with ",
                length(mz_values),
                " m/z values. The uploaded file will be used directly for m/z matching."
            )
        }, silent = TRUE)
    })

    shiny::observeEvent(input$toggle_search_sidebar, {
        bslib::toggle_sidebar(id = "search_sidebar", session = session)
    })

    shiny::observeEvent(input$run_search, {
        state$search_status <- "Running search..."
        requested_mode <- input$search_mode %||% "mass_match"
        search_grouping_context <- list(
            mode = requested_mode,
            uploaded = FALSE,
            query_count = NA_integer_
        )

        search_result <- tryCatch(
            {
                shiny::withProgress(message = "Running Search module", value = 0.1, {
                    adducts <- parse_adducts(
                        input$search_adducts,
                        default = default_search_adduct(input$search_ion_mode %||% "positive"),
                        allow_null = FALSE
                    )
                    mode <- requested_mode

                    if (mode == "target_mass_search") {
                        shiny::incProgress(0.3, detail = "Preparing target masses")
                        target_masses <- parse_numeric_values(input$search_target_masses)
                        if (is.null(input$search_target_feature_file)) {
                            stop("Please upload a feature table (first column mz, second column time) for target monoisotopic mass search.")
                        }
                        feature_table <- coerce_mz_time_columns(read_uploaded_table(input$search_target_feature_file))
                        shiny::incProgress(0.4, detail = "Matching target masses")
                        out <- mm_call(
                            "target_mass_search",
                            mass = target_masses,
                            feature_table = feature_table,
                            adduct = adducts,
                            mz_ppm = input$search_ppm
                        )
                    } else if (mode == "mass_match") {
                        shiny::incProgress(0.3, detail = "Preparing unknown features")
                        unknown_mz <- parse_text_or_upload_mz(
                            text_value = input$search_mz_values,
                            uploaded_file = input$search_unknown_file,
                            context_label = "search upload"
                        )
                        search_grouping_context$uploaded <- !is.null(input$search_unknown_file)
                        search_grouping_context$query_count <- length(unknown_mz)
                        unknown_feature <- mz_values_to_unknown_feature(unknown_mz)
                        prefiltered_database <- NULL
                        if (should_use_search_broad_prefilter(
                            mz_values = unknown_mz,
                            adducts = adducts,
                            uploaded_file = input$search_unknown_file
                        )) {
                            prefilter_detail <- if (!is.null(input$search_unknown_file)) {
                                "Broad-prefiltering database for uploaded m/z list"
                            } else {
                                "Broad-prefiltering database for multi-m/z/adduct search"
                            }
                            shiny::incProgress(0.15, detail = prefilter_detail)
                            prefiltered_database <- mm_call(
                                "prefilter_metabolite_database_by_mz_range",
                                mz_values = unknown_mz,
                                adducts = adducts,
                                ppm_threshold = input$search_ppm,
                                database = input$search_database
                            )
                        }
                        shiny::incProgress(0.4, detail = "Running mass_match")
                        out <- mm_call(
                            "mass_match",
                            unknown_feature = unknown_feature,
                            mz_ppm = input$search_ppm,
                            adduct = adducts,
                            database = input$search_database,
                            metabolite_database = prefiltered_database
                        )
                    } else {
                        shiny::incProgress(0.3, detail = "Preparing unknown features")
                        unknown_mz <- parse_text_or_upload_mz(
                            text_value = input$search_mz_values,
                            uploaded_file = input$search_unknown_file,
                            context_label = "search upload"
                        )
                        search_grouping_context$uploaded <- !is.null(input$search_unknown_file)
                        search_grouping_context$query_count <- length(unknown_mz)
                        unknown_feature <- mz_values_to_unknown_feature(unknown_mz)
                        prefiltered_database <- NULL
                        if (should_use_search_broad_prefilter(
                            mz_values = unknown_mz,
                            adducts = adducts,
                            uploaded_file = input$search_unknown_file
                        )) {
                            prefilter_detail <- if (!is.null(input$search_unknown_file)) {
                                "Broad-prefiltering database for uploaded m/z list"
                            } else {
                                "Broad-prefiltering database for multi-m/z/adduct search"
                            }
                            shiny::incProgress(0.15, detail = prefilter_detail)
                            prefiltered_database <- mm_call(
                                "prefilter_metabolite_database_by_mz_range",
                                mz_values = unknown_mz,
                                adducts = adducts,
                                ppm_threshold = input$search_ppm,
                                database = input$search_database
                            )
                        }
                        shiny::incProgress(0.4, detail = "Running isotopologue_mass_match")
                        out <- mm_call(
                            "isotopologue_mass_match",
                            unknown_feature = unknown_feature,
                            mz_ppm = input$search_ppm,
                            adduct = adducts,
                            database = input$search_database,
                            metabolite_database = prefiltered_database
                        )
                    }
                    shiny::incProgress(0.2, detail = "Done")
                    out
                })
            },
            error = function(e) {
                shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
                NULL
            }
        )

        if (!is.null(search_result)) {
            state$search_result <- search_result
            state$search_mode_used <- requested_mode
            state$search_grouping <- search_grouping_context
            state$search_run_stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
            state$search_status <- paste0("Search completed. Rows: ", nrow(search_result))
            session$onFlushed(function() {
                bslib::toggle_sidebar(id = "search_sidebar", open = FALSE, session = session)
            }, once = TRUE)
        } else {
            state$search_result <- NULL
            state$search_grouping <- NULL
            state$search_status <- "Search failed."
        }
    })

    output$search_table_single <- DT::renderDT({
        shiny::req(state$search_result)
        make_datatable(prepare_export_table(state$search_result))
    })

    output$search_group_page_size_ui <- shiny::renderUI({
        specs <- search_group_specs()

        if (!identical(specs$type, "grouped")) {
            return(NULL)
        }

        shiny::tags$div(
            style = "max-width:220px;margin-bottom:10px;",
            shiny::selectInput(
                inputId = "search_group_page_size",
                label = "Rows per group page",
                choices = c("10", "25", "50", "100"),
                selected = "10"
            )
        )
    })

    output$search_grouped_ui <- shiny::renderUI({
        specs <- search_group_specs()

        if (identical(specs$type, "empty")) {
            return(shiny::tags$em("No rows returned."))
        }

        if (identical(specs$type, "single")) {
            return(DT::DTOutput("search_table_single"))
        }

        page_size <- search_group_page_size()
        cards <- lapply(specs$groups, function(spec) {
            shiny::tags$div(
                class = "search-group-card",
                shiny::tags$div(
                    class = "search-group-header",
                    shiny::tags$strong(spec$heading)
                ),
                shiny::tags$div(
                    class = "search-group-body",
                    make_grouped_search_datatable(spec$data, page_size = page_size)
                )
            )
        })

        shiny::tagList(cards)
    })

    output$download_search <- shiny::downloadHandler(
        filename = function() {
            mode_label <- switch(
                state$search_mode_used %||% "search",
                mass_match = "mz-matching",
                isotopologue_mass_match = "isotopologue-mz-matching",
                target_mass_search = "target-monoisotopic-mass-search",
                "search"
            )
            stamp <- state$search_run_stamp %||% format(Sys.time(), "%Y%m%d-%H%M%S")
            paste0("massmatcher_", mode_label, "_", stamp, ".csv")
        },
        content = function(file) {
            shiny::req(state$search_result)
            utils::write.csv(prepare_export_table(state$search_result), file, row.names = FALSE)
        }
    )

    observe_cluster_file <- function() {
        if (isTRUE(state$cluster_use_demo)) {
            if (is.null(state$cluster_demo_table)) {
                state$cluster_demo_table <- load_demo_feature_table(rows = 500, cols = 12)
            }
            return(state$cluster_demo_table)
        }

        if (!is.null(input$cluster_feature_file)) {
            return(coerce_mz_time_columns(read_uploaded_table(input$cluster_feature_file)))
        }

        stop("Upload a feature table or click 'Use bundled demo features'.")
    }

    shiny::observeEvent(input$cluster_use_demo, {
        demo <- tryCatch(
            load_demo_feature_table(rows = 500, cols = 12),
            error = function(e) {
                shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
                NULL
            }
        )

        if (is.null(demo)) {
            return()
        }

        state$cluster_use_demo <- TRUE
        state$cluster_demo_table <- demo
        state$cluster_status <- "Bundled demo feature table selected for clustering."
    })

    shiny::observeEvent(input$cluster_feature_file, {
        shiny::req(input$cluster_feature_file)
        state$cluster_use_demo <- FALSE
        state$cluster_demo_table <- NULL
        state$cluster_status <- paste0("Uploaded feature table selected: ", input$cluster_feature_file$name)
    })

    shiny::observeEvent(input$run_cluster, {
        state$cluster_status <- "Running clustering..."
        cluster_error <- NULL

        clustering_result <- tryCatch(
            {
                shiny::withProgress(message = "Running Clustering module", value = 0, {
                    feature_table <- observe_cluster_file()
                    chosen_adducts <- parse_adducts(
                        input$cluster_adducts,
                        default = default_cluster_adducts(input$cluster_ion_mode %||% "positive"),
                        allow_null = TRUE
                    )

                    imputation <- if (identical(input$cluster_imputation, "None")) {
                        NA_character_
                    } else {
                        input$cluster_imputation
                    }

                    if (identical(input$cluster_database, "pubchem")) {
                        out <- mm_call(
                            "mz_match_clustering_pubchem",
                            met_raw_wide = feature_table,
                            mz_threshold = input$cluster_mz_threshold,
                            All_Adduct = chosen_adducts,
                            ion_mode = input$cluster_ion_mode,
                            adduct_correlation_r_threshold = input$cluster_adduct_r,
                            adduct_correlation_time_threshold = input$cluster_adduct_time,
                            isotopic_correlation_r_threshold = input$cluster_isotope_r,
                            isotopic_correlation_time_threshold = input$cluster_isotope_time,
                            imputation_method = imputation,
                            write_output = FALSE,
                            show_progress = FALSE,
                            progress_callback = function(step, total, detail) {
                                shiny::setProgress(value = step / total, detail = detail)
                            }
                        )
                    } else {
                        out <- mm_call(
                            "mz_match_clustering",
                            met_raw_wide = feature_table,
                            database = input$cluster_database,
                            mz_threshold = input$cluster_mz_threshold,
                            All_Adduct = chosen_adducts,
                            ion_mode = input$cluster_ion_mode,
                            adduct_correlation_r_threshold = input$cluster_adduct_r,
                            adduct_correlation_time_threshold = input$cluster_adduct_time,
                            isotopic_correlation_r_threshold = input$cluster_isotope_r,
                            isotopic_correlation_time_threshold = input$cluster_isotope_time,
                            imputation_method = imputation,
                            write_output = FALSE,
                            show_progress = FALSE,
                            progress_callback = function(step, total, detail) {
                                shiny::setProgress(value = step / total, detail = detail)
                            }
                        )
                    }
                    shiny::setProgress(value = 1, detail = "Done")
                    out
                })
            },
            error = function(e) {
                cluster_error <<- conditionMessage(e)
                shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
                NULL
            }
        )

        if (!is.null(clustering_result)) {
            state$clustering_result <- clustering_result
            table_names <- names(clustering_result)
            shiny::updateSelectInput(
                session = session,
                inputId = "cluster_view_table",
                choices = table_names,
                selected = table_names[[1]]
            )
            state$cluster_status <- "Clustering completed."
        } else {
            state$cluster_status <- paste0("Clustering failed: ", cluster_error %||% "unknown error")
        }
    })

    output$cluster_table <- DT::renderDT({
        shiny::req(state$clustering_result, input$cluster_view_table)
        make_datatable(state$clustering_result[[input$cluster_view_table]])
    })

    output$download_cluster <- shiny::downloadHandler(
        filename = function() {
            paste0("massmatcher_", input$cluster_view_table %||% "clustering", "_", Sys.Date(), ".csv")
        },
        content = function(file) {
            shiny::req(state$clustering_result)
            nm <- input$cluster_view_table %||% names(state$clustering_result)[[1]]
            if (!nm %in% names(state$clustering_result)) {
                nm <- names(state$clustering_result)[[1]]
            }
            utils::write.csv(as.data.frame(state$clustering_result[[nm]]), file, row.names = FALSE)
        }
    )

    output$download_cluster_all <- shiny::downloadHandler(
        filename = function() {
            paste0("massmatcher_clustering_all_", Sys.Date(), ".zip")
        },
        content = function(file) {
            shiny::req(state$clustering_result)

            temp_dir <- tempfile("massmatcher_cluster_export_")
            dir.create(temp_dir, recursive = TRUE, showWarnings = FALSE)
            on.exit(unlink(temp_dir, recursive = TRUE, force = TRUE), add = TRUE)

            csv_files <- character(0)
            for (nm in names(state$clustering_result)) {
                path <- file.path(temp_dir, paste0(nm, ".csv"))
                utils::write.csv(as.data.frame(state$clustering_result[[nm]]), path, row.names = FALSE)
                csv_files <- c(csv_files, path)
            }

            old_wd <- getwd()
            on.exit(setwd(old_wd), add = TRUE)
            setwd(temp_dir)
            utils::zip(zipfile = file, files = basename(csv_files))
        }
    )

    output$download_cluster_demo <- shiny::downloadHandler(
        filename = function() {
            paste0("feature_table_exp_hilicpos_demo_", Sys.Date(), ".csv")
        },
        content = function(file) {
            demo <- state$cluster_demo_table
            if (is.null(demo)) {
                demo <- load_demo_feature_table(rows = 500, cols = 12)
            }
            utils::write.csv(as.data.frame(demo), file, row.names = FALSE)
        }
    )

    shiny::observeEvent(input$run_metadata, {
        state$metadata_status <- "Preparing Additional Meta Data enrichment..."
        update_metadata_progress(
            value = 2,
            detail = "Preparing Additional Meta Data enrichment",
            animated = TRUE
        )

        metadata_setup <- tryCatch(
            {
                shiny::req(input$meta_upload_file)
                uploaded <- read_uploaded_table(input$meta_upload_file)
                list(
                    uploaded = uploaded,
                    database = input$meta_database,
                    biospecimen = input$meta_biospecimen
                )
            },
            error = function(e) {
                shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
                NULL
            }
        )

        if (is.null(metadata_setup)) {
            state$metadata_status <- "Metadata enrichment failed."
            update_metadata_progress(
                value = 0,
                detail = "Metadata enrichment failed to start.",
                label = "Additional Meta Data progress",
                striped = FALSE,
                animated = FALSE
            )
            return()
        }

        metadata_result <- tryCatch(
            {
                progress <- shiny::Progress$new(session, min = 0, max = 1)
                on.exit(progress$close(), add = TRUE)

                progress$set(
                    message = "Running Additional Meta Data module",
                    value = 0.05,
                    detail = "Reading uploaded table"
                )
                state$metadata_status <- "Reading uploaded table..."
                update_metadata_progress(
                    value = 5,
                    detail = "Reading uploaded table",
                    animated = TRUE
                )

                out <- mm_call(
                    "annotate_match_results",
                    result_table = metadata_setup$uploaded,
                    database = metadata_setup$database,
                    biospecimen = metadata_setup$biospecimen,
                    query_missing_classification = TRUE,
                    progress_callback = function(current, total, detail) {
                        fraction <- if (!is.finite(total) || total <= 0) {
                            0
                        } else {
                            min(max(current / total, 0), 1)
                        }
                        progress$set(
                            message = "Running Additional Meta Data module",
                            value = 0.15 + (0.75 * fraction),
                            detail = detail
                        )
                        update_metadata_progress(
                            value = round((0.15 + (0.75 * fraction)) * 100, 0),
                            detail = detail,
                            animated = TRUE
                        )
                    },
                    status_callback = function(detail) {
                        state$metadata_status <- detail
                    }
                )

                state$metadata_status <- "Finalizing enriched metadata table..."
                progress$set(
                    message = "Running Additional Meta Data module",
                    value = 0.95,
                    detail = "Finalizing enriched metadata table"
                )
                update_metadata_progress(
                    value = 95,
                    detail = "Finalizing enriched metadata table",
                    animated = TRUE
                )
                out
            },
            error = function(e) {
                update_metadata_progress(
                    value = 100,
                    detail = conditionMessage(e),
                    label = "Additional Meta Data failed",
                    striped = FALSE,
                    animated = FALSE
                )
                shiny::showNotification(conditionMessage(e), type = "error", duration = 10)
                NULL
            }
        )

        if (is.null(metadata_result)) {
            state$metadata_status <- "Metadata enrichment failed."
            return()
        }

        state$metadata_table <- metadata_result
        state$metadata_status <- paste0("Metadata enrichment completed. Rows: ", nrow(metadata_result))
        update_metadata_progress(
            value = 100,
            detail = paste0("Metadata enrichment completed. Rows: ", nrow(metadata_result)),
            label = "Additional Meta Data complete",
            striped = FALSE,
            animated = FALSE
        )
    })

    output$metadata_table <- DT::renderDT({
        shiny::req(state$metadata_table)
        make_datatable(state$metadata_table)
    })

    output$download_metadata <- shiny::downloadHandler(
        filename = function() {
            paste0("massmatcher_metadata_", Sys.Date(), ".csv")
        },
        content = function(file) {
            shiny::req(state$metadata_table)
            utils::write.csv(prepare_export_table(state$metadata_table), file, row.names = FALSE)
        }
    )
}

shiny::shinyApp(ui = ui, server = server)
