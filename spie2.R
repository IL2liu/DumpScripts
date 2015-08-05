addspie <- function(expected, observed, x, y, R, ...) {
# Plot a Spie chart

   # A bit of maths.
   angle <- 2*pi * diffinv(expected / sum(expected))
   area <- pi * observed / sum(observed)
   radius <- R * sqrt(2 * area / diff(angle))

   # Manually recycle the parameters passed to 'polygon()'.
   # NB: 'polargs.i' is a list of lists of arguments.
   dotargs = list(...)
   polargs <- rep(
      do.call(
         what = mapply,
         args = c(list(FUN=list, SIMPLIFY=FALSE), dotargs)
      ),
      length.out = length(radius)
   )

   for (i in 1:length(radius)) {
      nsteps <- (angle[i+1] - angle[i]) / .01
      steps <- seq(from=angle[i], to=angle[i+1], length.out=nsteps)
      x_ <- c(x, x + radius[i]*cos(steps))
      y_ <- c(y, y + radius[i]*sin(steps))
      do.call(what=polygon, args=c(list(x=x_, y=y_), polargs[[i]]))
   }
   steps <- seq(0, 2*pi, length.out=628)
   lines(x + R* cos(steps), y + R * sin(steps), lty=2)
}
