# # Image Classification With Cubical Filtrations and Persistence Images

# In this example, we will show how to use Ripserer in an image classification
# context. Persistent homology is not a predictive algorithm, but it can be used to extract
# useful features from data.

using Ripserer
using PersistenceDiagrams
using Images # also required: ImageIO to read .png files
using Plots
using ProgressMeter
using Random
Random.seed!(1337)

data_dir = joinpath(@__DIR__, "../assets/data/malaria") # replace with the correct path.
nothing # hide

# Let's load the data. We will use a a [data
# set](https://lhncbc.nlm.nih.gov/publication/pub9932) with microscope images of healthy
# cells and cells infected with malaria. The original data set is quite large, but we can
# pretend we were only given 200 images to work with. We have chosen the 200 images
# randomly.

uninfected = shuffle!(load.(readdir(joinpath(data_dir, "uninfected"), join=true)))
infected = shuffle!(load.(readdir(joinpath(data_dir, "infected"), join=true)))

images = [uninfected; infected]
isinfected = [fill(false, length(uninfected)); fill(true, length(infected))]
nothing # hide

# Let's see what the images look like.

plot(plot(uninfected[1], title="Healthy cell"),
     plot(uninfected[2], title="Healthy cell"),
     plot(infected[1], title="Infected cell"),
     plot(infected[2], title="Infected cell"))

# Next, we convert images to floating point gray scale values. We do not have to resize the
# images. Maybe some additional preprocessing, such as normalization would help, but we'll
# skip it for this example.

inputs = [Float32.(Gray.(image)) for image in images]
nothing # hide

# Now we can compute persistence diagrams. Since we are working with images, we have to use
# the `Cubical` filtration type. Cubical persistent homology should detect the dark spots
# (local minima) in the images. It's pretty efficient, so this should only take a few
# seconds.

diagrams = @showprogress [ripserer(Cubical(i)) for i in inputs]
nothing # hide

# This is what some of the diagrams look like.

plot(plot(images[1], title="Healthy"), plot(diagrams[1]))

#

plot(plot(images[end], title="Infected"), plot(diagrams[end]))

# Notice that there is a lot more going on in the middle of the infected diagram, especially
# in ``H_0``.

# The persistence diagrams are nice, but are hard to use with machine learning algorithms
# since the number of points in the diagram may be different for every image. We can solve
# this problem by using a vectorization method, such as converting all diagrams to
# persistence images.

# Persistence images work by weighting each point in the diagram with a distribution. The
# distribution defaults to a Gaussian, but any function of two arguments can be used. Each
# point is also weighted by a weighting function that should be equal to zero along the
# ``x``-axis. It defaults to a function that is zero on the ``x``-axis and linearly
# increases to the maximum persistence in the diagram.

# We start by splitting the diagrams into their 0 and 1 dimensional components.

dim_0 = first.(diagrams)
dim_1 = last.(diagrams)
nothing; # hide

# We feed the diagram to the `PersistenceImage` constructor which will choose ranges that
# will fit all the diagrams. We set the sigma value to 0.1, since all persistence pairs are
# in the ``[0,1]×[0,1]`` square and the default sigma of 1 would be too wide.

image_0 = PersistenceImage(dim_0, size=(5, 5), sigma=0.1)

#

image_1 = PersistenceImage(dim_1, size=(5, 5), sigma=0.1)

#

plot(plot(dim_0[end], persistence=true),
     heatmap(image_0(dim_0[end]), aspect_ratio=1))

#

plot(plot(dim_1[end], persistence=true),
     heatmap(image_1(dim_1[end]), aspect_ratio=1))

# We convert the diagrams to images and use `vec` to turn them into flat vectors. We then
# concatenate them. The result is a vector of length 50.

persims = [[vec(image_0(dim_0[i])); vec(image_1(dim_1[i]))]
           for i in 1:length(diagrams)]

# Now it's time to fit our model. We will use
# [GLMNet.jl](https://github.com/JuliaStats/GLMNet.jl) to fit a regularized linear model.

using GLMNet

# Convert the image vectors to a matrix that will be understood by `glmnet`.

X = reduce(hcat, persims)'
y = outputs
nothing; # hide

# Start by randomly splitting the data into two sets, a training and a testing set.

perm = shuffle(1:200)
train_x = X[perm[1:100], :]
train_y = y[perm[1:100]]
test_x = X[perm[101:end], :]
test_y = y[perm[101:end]]
nothing; # hide

# Fit the model and predict.

path = glmnet(train_x, train_y)
cv = glmnetcv(train_x, train_y)

λ = path.lambda[argmin(cv.meanloss)]
path = glmnet(train_x, train_y; lambda=[λ])

predictions = .!iszero.(round.(GLMNet.predict(path, test_x)))
nothing; # hide

# Get the classification accuracy.

accuracy = count(predictions .== test_y) / length(test_y)

# Not half bad considering we haven't touched the images and we left pretty much all
# settings on default.

# Now let's look at the misclassified examples.

missed = findall(predictions .!= test_y)
label = ("Healthy", "Infected")
plts = [plot(images[i],
             title="$(label[test_y[i] + 1])",
             ticks=nothing)
        for i in missed]
plot(plts...)

# Finally, let's look at which parts of the persistence images `glmnet` considered important.

plot(heatmap(reshape(path.betas[1:25], (5,5)), title="H₀ coefficients"),
     heatmap(reshape(path.betas[26:50], (5,5)), title="H₁ coefficients"))

# These correspond to the area we identified at the beginning. Also note that in this case,
# the classifier does not care about ``H_1`` at all.
