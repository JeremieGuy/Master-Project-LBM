import imageio


directory = "Monitoring/Simple_tube_square_obstacle_5000_it/FoutEvo_node[130,50]"
maxIter = 5000
plots = 10

frames = []
t = maxIter//plots

print("Making Gif ...")
for i in range(t):
    # print(i)
    num = "{0:0=5d}".format(i)
    image = imageio.v2.imread(f"./" + directory + "/node_" + str(num) + ".png")
    frames.append(image)

imageio.mimsave("./" + directory + "/_nodeEvo.gif", frames, duration = 80)

print("\nDone.")