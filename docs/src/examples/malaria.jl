# # Image Classification With Cubical Persistent Homology

# In this example, we will show how to use Ripserer in an image classification
# context. Persistent homology is not a predictive algorithm, but it can be used to extract
# useful features from data.

# ## Setting up

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

uninfected = shuffle!(load.(readdir(joinpath(data_dir, "uninfected"); join=true)))
infected = shuffle!(load.(readdir(joinpath(data_dir, "infected"); join=true)))

images = [uninfected; infected]
classes = [fill(false, length(uninfected)); fill(true, length(infected))]
nothing # hide

# Let's see what the images look like.

plot(
    plot(uninfected[1]; title="Healthy"),
    plot(uninfected[2]; title="Healthy"),
    plot(infected[1]; title="Infected"),
    plot(infected[2]; title="Infected"),
)

# To make the images work with Ripserer, we convert them to floating gray scale values. We
# do not have to resize the images. Maybe some additional preprocessing, such as
# normalization would help, but we'll skip it for this example.

inputs = [Gray.(image) for image in images]
nothing # hide

# Now we can compute persistence diagrams. Since we are working with images, we have to use
# the `Cubical` filtration type. Cubical persistent homology should detect the dark spots
# (local minima) in the images. It's pretty efficient, so this should only take a few
# seconds.

diagrams = @showprogress [ripserer(Cubical(i)) for i in inputs]

# This is what some of the diagrams look like.

plot(plot(images[1]; title="Healthy"), plot(diagrams[1]))

#

plot(plot(images[end]; title="Infected"), plot(diagrams[end]))

# Notice that there is a lot more going on in the middle of the infected diagram, especially
# in ``H_0``.

# The persistence diagrams might look nice, but are hard to use with machine learning
# algorithms. The number of points in the diagram may be different for every image, even
# when images are of the same size. We can solve this problem by using a vectorization
# method, such as converting all diagrams to persistence images.

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
# in the ``[0,1]×[0,1]`` square and the default sigma of 1 would be too wide. We will use
# the default image size, which is 5×5.

image_0 = PersistenceImage(dim_0; sigma=0.1)

#

image_1 = PersistenceImage(dim_1; sigma=0.1)

# Let's see how some of the images look like.

plot(plot(dim_0[end]; persistence=true), heatmap(image_0(dim_0[end]); aspect_ratio=1))

#

plot(plot(dim_1[end]; persistence=true), heatmap(image_1(dim_1[end]); aspect_ratio=1))

# Next, we convert all diagrams to images and use `vec` to turn them into flat vectors. We
# then concatenate the zero and one-dimensional images. The result is a vector of length 50
# for each diagram.

persims = [[vec(image_0(dim_0[i])); vec(image_1(dim_1[i]))] for i in 1:length(diagrams)]

# ## Fitting A Model

# Now it's time to fit our model. We will use
# [GLMNet.jl](https://github.com/JuliaStats/GLMNet.jl) to fit a regularized linear model.

using GLMNet

# Convert the image vectors to a matrix that will be understood by `glmnet`.

X = reduce(hcat, persims)'
y = classes
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

count(predictions .== test_y) / length(test_y)

# Not half bad considering we haven't touched the images and we left pretty much all
# settings on default.

# Now let's look at the misclassified examples.

missed = findall(predictions .!= test_y)
label = ("Healthy", "Infected")
plts = [plot(images[i]; title="$(label[test_y[i] + 1])", ticks=nothing) for i in missed]
plot(plts...)

# Finally, let's look at which parts of the persistence images `glmnet` considered important.

plot(
    heatmap(reshape(path.betas[1:25], (5, 5)); title="H₀ coefficients"),
    heatmap(reshape(path.betas[26:50], (5, 5)); title="H₁ coefficients"),
)

# These correspond to the area we identified at the beginning. Also note that in this case,
# the classifier does not care about ``H_1`` at all.

# ## Using MLJ

# Another, more straightforward way to execute a similar pipeline is to use Ripserer's
# [MLJ.jl](https://github.com/alan-turing-institute/MLJ.jl) integration. We will use a
# random forest classifier for this example.

# We start by loading MLJ and the classifier. Not that the
# [DecisionTree.jl](https://github.com/bensadeghi/DecisionTree.jl) package needs to be
# installed for this to work.

using MLJ
tree = @load RandomForestClassifier pkg = "DecisionTree" verbosity = 0

# We create a pipeline of `CubicalPersistentHomology` followed by the classifier. In this
# case, `CubicalPersistentHomology` takes care of both the homology computation and the
# conversion to persistence images.

pipe = @pipeline(CubicalPersistentHomology(), tree)

# We train the pipeline the same way you would fit any other MLJ model. Remember, we need to
# use grayscale versions of images stored in `inputs`.

classes = coerce(classes, Binary)
train, test = partition(eachindex(classes), 0.7; shuffle=true, rng=1337)
mach = machine(pipe, inputs, classes)
fit!(mach; rows=train)

# Next, we predict the classes on the test data and print out the classification accuracy.

yhat = predict_mode(mach, inputs[test])
accuracy(yhat, classes[test])

# The result is quite a bit worse than before. We can try mitigating that by using a
# different vectorizer.

pipe.cubical_persistent_homology.vectorizer = PersistenceCurveVectorizer()
mach = machine(pipe, inputs, classes)
fit!(mach; rows=train)

yhat = predict_mode(mach, inputs[test])
accuracy(yhat, classes[test])

# The result could be improved further by choosing a different model and
# vectorizer. However, this is just a short introduction. Please see the [MLJ.jl
# documentation](https://alan-turing-institute.github.io/MLJ.jl/dev/) for more information
# on model tuning and selection, and the [PersistenceDiagrams.jl
# documentation](https://mtsch.github.io/PersistenceDiagrams.jl/dev/mlj/) for a list of
# vectorizers and their options.
