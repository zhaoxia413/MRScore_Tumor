library(remotes)
devtools::install_github("EmilHvitfeldt/miscpalettes")
library(miscpalettes)
 pals::pal.bands(
  artistic[[1]],
  artistic[[2]],
  artistic[[3]],
  artistic[[4]],
  artistic[[5]],
  labels = names(artistic),
  main = "artistic palettes"
)
pals::pal.bands(
  mschart[[1]],
  mschart[[2]],
  mschart[[3]],
  mschart[[4]],
  mschart[[5]],
  mschart[[6]],
  mschart[[7]],
  mschart[[8]],
  mschart[[9]],
  mschart[[10]],
  mschart[[11]],
  mschart[[12]],
  labels = names(mschart),
  main = "mschart palettes"
)
library(devtools)  
library(lisa)
par(mfrow = c(6, 3))
lapply(sample(lisa, 18), plot)

x <- lisa_palette("JackBush_1", 1000, "continuous")
y <- lisa_palette("PabloPicasso", 2, "discrete")
z <- lisa_palette("KatsushikaHokusai", 1000, "continuous")
lapply(list(x, y, z), plot)
useage: scale_color_manual(values = lisa$`Jean-MichelBasquiat`)

devtools::install_github("jaredhuling/jcolors")
library(jcolors)
jcolors('default')
display_all_jcolors()
display_all_jcolors_contin()
display_jcolors("pal4")
display_jcolors("pal12")
grid.arrange(pltl + scale_color_jcolors(palette = "default"),
             pltd + scale_color_jcolors(palette = "default"), ncol = 2)
install.packages("harrypotter")
library(harrypotter)
pal <- hp(n = 8, house = "Gryffindor")
image(volcano, col = pal)
pal <- hp(n = 128, house = "Gryffindor")
image(volcano, col = pal)
scale_fill_hp(house = "hufflepuff")
scale_fill_hp(discrete = TRUE, option = "ronweasley2", name = "Cut")
scale_colour_hp_d(option = "LunaLovegood", name = "Clarity")
install.packages('ghibli')
library(ghibli)
par(mfrow=c(9,3))
for(i in names(ghibli_palettes)) print(ghibli_palette(i))
ghibli_palettes$MarnieLight1
ghibli_palettes$PonyoLight
scale_colour_ghibli_d("LaputaMedium", direction = -1)

install.packages("gameofthrones")
library(gameofthrones)
pal <- got(250, option = "Targaryen")
image(volcano, col = pal)
scale_fill_got(option = "Martell")
scale_fill_got(option = "Baratheon")
scale_fill_got_d(option = "Daenerys", direction = - 1) 
scale_fill_got(discrete = TRUE, option = "Margaery")

install.packages("fishualize")
library(fishualize)
fishualize::fish_palettes()
scale_fill_fish_d(option = "Scarus_quoyi")
scale_color_fish(option = "Lepomis_megalotis", direction = -1)
scale_fill_fish(option = "Zebrasoma_velifer", trans = "sqrt")


library(DresdenColor)
devtools::install_github("katiesaund/DresdenColor")
dresden_palette("briefcases", type = "discrete")
dresden_palette("briefcases", n = 50, type = "continuous")
dresden_palette(type = "paired")
mat <- matrix(rnorm(n = 100, mean = 0, sd = 20), ncol = 10)
heatmap(mat, col = dresden_palette("foolmoon", type = "continuous", n = length(unique(mat))))

devtools::install_github("an-bui/calecopal")
library(calecopal)
# all palettes
names(cal_palettes)
cal_palette(name = "desert", n = 15, type = "continuous")
cal_palette("sierra1", n = 50, type = "continuous")
scale_fill_manual(values = cal_palette("sierra1"))

install.packages("basetheme")
library(basetheme)
basetheme(pch=19, mgp=c(2,.7,0), tck=-.01)
basetheme("clean")
basetheme("brutal")
plot(hclust(dist(USArrests), "ward.D2"), hang=-1)

devtools::install_github('awhstin/awtools')
devtools::install_github("jakelawlor/PNWColors") 
library(PNWColors)
names(pnw_palettes)
#[1] "Starfish" "Shuksan"  "Bay"      "Winter"   "Lake"     "Sunset"   "Shuksan2" 
#[8] "Cascades" "Sailboat" "Moth" "Spring"   "Mushroom" "Sunset2"  "Anemone"    
pnw_palette(name="Starfish",n=7,type="discrete")
pnw_palette("Winter",100)
pnw_palette("Bay",8,type="continuous")
pnw_palette("Moth",12)
pal=pnw_palette("Lake",5, type = "discrete")
scale_fill_gradientn(colours = pal)
scale_fill_gradientn(colours = rev(pal))

devtools::install_github("m-clark/NineteenEightyR")
library(NineteenEightyR)

install.packages("yarrr")
library(yarrr)
yarrr::yarrr.guide()

pirateplot(formula = weight ~ Time,
           data = ChickWeight,
           main = "Chicken weights by Time (week)")
piratepal("all")
my.cols <- piratepal(palette = "google",  
                     trans = .5)
install.packages("palettesForR")
library(palettesForR)
names(palettesForR)
palettesForR::Android_gpl
palettesForR::Blues_gpl
palettesForR::Inkscape_gpl[1:10]
palettesForR::showPalette(Inkscape_gpl)
#很有用但是安不上,多试几次就按上了
devtools::install_github("m-clark/visibly")
install.packages("colortools")
library(colortools)
library(visibly)
library(plotly)
mycol=create_palette('papayawhip')
mycol
theme_plotly()
data('bfi', package = 'visibly')
cor_matrix = cor(bfi, use='pair')
corr_heat(cor_matrix)
fit_lm = lm(mpg ~ ., mtcars)
plot_coefficients(fit_lm)
library(mgcv)#
d = gamSim()
#Gu & Wahba 4 term additive model

data<-data%>%
  mutate(y=factor(Group)%>%
           fct_recode("2"="I_F",
                      "1"="I_N",
                      "0"="N",
                     "3"="P_N",
                      "4"="P_F"))
dfmat<-data[,-1]
dfvars=paste0(colnames(dfmat)[-21],collapse = "+")
fmla=as.formula(paste0("y~",dfvars))
gam_model=gam(fmla,data=dfmat)
gam_model = gam(y ~ x0 + s(x1) + s(x2, bs='gp') + s(x3, bs='ps'), data=d)
plot_gam(gam_model, main_var = Clostridium)
plot_gam(gam_model, main_var = x2)
plot_gam_check(gam_model)
fit_lm = lm(y ~ ., dfmat)
plot_coefficients(fit_lm)
install.packages("rcartocolor")
library(rcartocolor)
display_carto_all()
my_colors = carto_pal(7, "Burg")
my_colors
display_carto_all(colorblind_friendly = TRUE)
scale_fill_carto_c(name = "Life expectancy: ",
                   type = "diverging", palette = "Earth", direction = -1)

scale_fill_carto_d(name = "Region: ", palette = "Safe")

devtools::install_github("thomasp85/scico")
#heatmap
devtools::install_github("thomasp85/scico")
library(scico)
scico_palette_show()
scico(30, palette = 'lapaz')
scale_fill_scico(palette = 'davos') 

install.packages("trekcolors")
library(trekcolors)
trek_pal()
view_trek_pals()
trek_pal("starfleet")
p + scale_color_trek()
scale_color_lcars("2357") # equivalent to scale_color_trek("lcars_2357")
scale_fill_trek("klingon")
scale_fill_trek("romulan")
scale_fill_lcars("2357")
lcars_colors()
#Custom palettes
lcars_custom_pal <- lcars_colors_pal(c("pale-canary", "eggplant"))
lcars_custom_pal(8)
scale_fill_lcars2("pale-canary", "danub")
scale_fill_lcars2("pale-canary", "danub", dark = TRUE, divergent = TRUE)
install.packages("tvthemes")
library(tvthemes)
theme_brooklyn99(title.font = "Titillium Web",
                 text.font = "Calibri Light",
                 subtitle.size = 14)
scale_color_westeros(palette = "Targaryen")
## patchwork together:
arryn + manderly - martell + plot_layout(ncol = 1)




















