using SparseGridsKit, Luxor
Drawing(500, 500, "docs/src/assets.svg")
origin(250,250)
setcolor("red")
pts = get_grid_points(sg);
for i in 1:get_n_grid_points(sg)
    circle(Point(200*pts[i][1], 200*pts[i][2]), 20, action = :fill)
end
finish()