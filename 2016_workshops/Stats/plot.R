# Plot CCLE pharmacological data
plot.drug.activity <- function(
	pharm,
	compounds=NULL,
	samples=NULL,
	group=NULL,
	x.var="ec50",
	y.var="act.max",
	xlab=NULL,
	ylab=NULL,
	xlim=c(0.0025, 8),
	ylim=c(-100, 0),
	xat=c(0.0025, 0.008, 0.025, 0.08, 0.25, 0.8, 2.5, 8),
	pch=19,
	add=FALSE,
	log="x",
	legend.x="topleft",
	legend.y=NULL,
	col='blue',
	palette=c("#999999", "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF" ),
	draw.legend=TRUE,
	...
) {
	# set options
	#scipen <- 100;
	#scipen.old <- getOption("scipen");
	#options(scipen=scipen);

	# select samples and compounds
	if (is.null(samples)) {
		samples <- rep(TRUE, nrow(pharm[[x.var]]));
	}
	if (is.null(compounds)) {
		compounds <- rep(TRUE, ncol(pharm[[x.var]]));
	}
	x <- pharm[[x.var]][samples, compounds];
	y <- pharm[[y.var]][samples, compounds];

	# set axis labels
	if (is.null(xlab)) {
		xlab <- x.var;
	}
	if (is.null(ylab)) {
		ylab <- y.var;
	}

	# set colours
	if (!is.null(group)) {
		if (!is.factor(group)) {
			group <- factor(group);
		}
		n.groups <- length(unique(group));
		n.colours <- length(palette);
		if (n.groups > n.colours) {
			# too many groups: collapse groups into "other"
			# use the first colour for other
			palette <- palette[c(2:n.colours, 1)];
			group.int <- .bound(as.integer(group), c(1, n.colours));
			group <- factor(group.int, levels=1:n.colours, labels=c(levels(group)[1:(n.colours-1)], "other"));
			n.groups <- n.colours;
		}
		col.palette <- palette[1:n.groups];
		col <- .colour.by.factor(group, col=col.palette);
	}

	# plot
	if (add) {
		points(x, y, pch=pch, col=col, ...);
	} else {
		plot(x, y, xlim=xlim, ylim=ylim, pch=pch, xlab=xlab, ylab=ylab, col=col, log=log, xaxt="n", ...);
		axis(1, at=xat);
	}
	if (draw.legend && length(unique(col)) > 1 && !is.null(group)) {
		legend(legend.x, legend.y, levels(group), col = col.palette, bty="n", pch=pch);
	}

	# restore options
	#options(scipen=scipen.old);
}

# Represent each level of a factor as a colour
.colour.by.factor <- function(f, col) {
	as.character(factor(f, labels=col))
}

.bound <- function(x, lim) {
	x[x < lim[1]] <- lim[1];
	x[x > lim[2]] <- lim[2];
	x
}

