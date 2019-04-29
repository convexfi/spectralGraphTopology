# See https://www.zevross.com/blog/2017/06/19/tips-and-tricks-for-working-with-images-and-figures-in-r-markdown-documents/
# for info on resolution of figures and images

# To optimize png's, first install:
brew install optipng
brew install pngquant

```{r}
# Set up the hook
library(knitr)
knit_hooks$set(optipng = hook_optipng)
knit_hooks$set(pngquant = hook_pngquant)
```

Then use it in chunks:
```{r, optipng = "-o7", pngquant = "--speed=1"}
# Maximum optimization, size is 17kb
ggplot(cars, aes(speed, dist)) + geom_point()
```

#See this comparison: http://pointlessramblings.com/posts/pngquant_vs_pngcrush_vs_optipng_vs_pngnq/

# To compress a .png from command window:
optipng -o7 *.png
pngquant --speed=1 *.png  (this one is better)
