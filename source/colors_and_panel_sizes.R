H_FULL <- 11 / 4.5
H_HALF <- 11 / 9

W_FULL <- 8.5 / 3
W_DBL <- 8.5 / 3 * 2

library(paletteer)

df <- data.frame(
    palette_name = c("lisa::OskarSchlemmer", "NineteenEightyR::miami1", "LaCroixColoR::Lemon", "LaCroixColoR::Orange", "LaCroixColoR::CranRaspberry", "colorBlindness::Blue2Orange8Steps", "khroma::BuRd"),
    fill_NP = c(2, 2, 5, 5, 5, 2, 2),
    fill_RP = c(4, 4, 2, 2, 2, 7, 7),
    color_NP = c(1, 1, 6, 6, 6, 1, 1),
    color_RP = c(5, 5, 1, 1, 1, 8, 8),
    center = c(3, 3, NA, NA, NA, NA, NA)
)

RP_color_pallette <- "lisa::OskarSchlemmer"
NP_color_palette <- "LaCroixColoR::Lemon"
fill_NP <- paletteer_d(NP_color_palette)[[5]] |>
    as.character() |>
    substr(1, 7)
color_NP <- paletteer_d(NP_color_palette)[[6]] |>
    as.character() |>
    substr(1, 7)

fill_RP <- paletteer_d(RP_color_pallette)[[4]] |>
    as.character() |>
    substr(1, 7)
color_RP <- paletteer_d(RP_color_pallette)[[5]] |>
    as.character() |>
    substr(1, 7)

center <- "white"

color_AFR_high <- "#86325F"
color_AFR_low <- "#5A90C7"

fill_AFR_high <- "#9b4472"
fill_AFR_low <- "#699ed4"



fill_theme_binary <- c(
    "NP" = fill_NP,
    "RP" = fill_RP
)

color_theme_binary <- c(
    "NP" = color_NP,
    "RP" = color_RP
)

fill_theme_emphasis <- c(
    "NP" = fill_NP,
    "NP*" = color_NP,
    "RP" = fill_RP,
    "RP*" = color_RP,
    "No Difference" = "white"
)

fill_theme_binary_race <- c(
    "white" = fill_NP,
    "black_african_american" = fill_RP
)

fill_theme_binary_afr_binned <- c(
    "AFR Low" = fill_AFR_low,
    "AFR High" = fill_AFR_high
)

color_theme_binary_race <- c(
    "white" = color_NP,
    "black_african_american" = color_RP
)

color_theme_binary_afr_binned <- c(
    "AFR Low" = color_AFR_low,
    "AFR High" = color_AFR_high
)

fill_theme_emphasis_race <- c(
    "white" = fill_NP,
    "white*" = color_NP,
    "black_african_american" = fill_RP,
    "black_african_american*" = color_RP,
    "No Difference" = "white"
)

gradientn_colors <- c(fill_theme_binary[["RP"]], "white", fill_theme_binary[["NP"]])

fill_theme_continuous <- list(low = fill_NP, mid = "white", high = color_RP)

PAGE_WIDTH <- 210
PAGE_HEIGHT <- 297
PAGE_UNITS <- "mm"
