using MLJ
using Plots
using PersistenceDiagrams
using Ripserer

using Ripserer.MLJRipserer
using PersistenceDiagrams.MLJPersistenceDiagrams

# We define the linked twist map.

function linked_twist_map(n, r, x1=rand(), y1=rand())
    xs = zeros(n)
    ys = zeros(n)
    xs[1] = x1
    ys[1] = y1

    for i in 2:n
        xs[i] = mod(xs[i - 1] + r * ys[i - 1] * (1 - ys[i - 1]), 1)
        ys[i] = mod(ys[i - 1] + r * xs[i] * (1 - xs[i]), 1)
    end

    return unique!(collect(zip(xs, ys)))
end

# Next, we generate the data set. We will use the following parameters.

parameters = (2.5, 3.5, 4.0, 4.1, 4.3);

# For each parameter, we will generate 50 point clouds with 1000 points each. Note that
# since we're doing classification, we have to make the y values categorical. We will split
# the data into a tran set and a test set using the standard 70/30 split.

data = [(r, linked_twist_map(1000, r)) for r in parameters for _ in 1:50]
x = last.(data)
y = categorical(string.(first.(data)))

train, test, validate = partition(eachindex(y), 0.6, 0.2; shuffle=true, rng=1337);

# Let's look at a few examples.

plots = [
    scatter(x[i]; title="r=$(y[i])", markersize=0, legend=false, aspect_ratio=1)
    for i in [1, 51, 101, 151, 201, 2, 52, 102, 152, 202]
]
plot(plots...; layout=(2, 5), size=(1000, 400))

# Now, we can create a machine learning pipeline. Ripserer and PersistenceDiagrams provide
# two types of unsupervised models.

# The first is used to compute persistent homology. This
# currently includes `RipsPersistentHomology`, `AlphaPersistentHomology`, and
# `CubicalPersistentHomology` which correspond to to `Rips`, `Alpha`, and `Cubical`
# filtrations defined by `Ripserer`.

# Since the output of persistent homology is a collection of `PersistenceDiagram`s, we need
# another set of models that can translate those into vectors of numbers, which can then be
# used as features by other models in the MLJ ecosystem. The list currently includes
# `PersistenceImageVectorizer`, `PersistenceCurveVectorizer`, and
# `PersistenceLandscapeVectorizer`. Please see [the PersistenceDiagrams.jl
# documentation](TODO) for the descriptions of their hyperparameters.

# In this example, we will focus on persistent homology, so we will use a simple random
# forest classifier to make predictions.

tree = @load RandomForestClassifier pkg = "DecisionTree"

pipe = @pipeline(
    AlphaPersistentHomology(),
    PersistenceImageVectorizer(; width=10, height=10, sigma=0.05),
    tree,
)

mach = machine(pipe, x, y)
@time fit!(mach; rows=train);
@profview yhat = predict_mode(mach, x[test]);
accuracy(yhat, y[test])

curve_pipe = @pipeline(AlphaPersistentHomology(), PersistenceCurveVectorizer(), tree,)

curves = [:betti_curve, :life, :midlife, :life_entropy, :midlife_entropy]
curve_range = range(curve_pipe, :(persistence_curve_vectorizer.curve); values=curves)
length_range = range(
    curve_pipe, :(persistence_curve_vectorizer.length); values=[3, 5, 10, 15, 20, 40, 100]
)

tuned_curves = TunedModel(;
    model=curve_pipe, ranges=[curve_range, length_range], measure=cross_entropy
)
#tuned_curves.tuning = Grid(goal=100)
mach = machine(tuned_curves, x, y)
fit!(mach; rows=train)
r = report(mach)
r.best_model

res = r.plotting
vals_sf = res.parameter_values[:, 1]
vals_bf = res.parameter_values[:, 2]

yhat = predict_mode(mach, x[test])

accuracy(yhat, y[test])
