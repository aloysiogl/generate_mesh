import mesh

# mesh.generateMesh([1, 1, 1], [0.1, 0.5, 0.5], [0.1, 0.5, 0.5], [0.2, 0.1, 0.25], [1, 1, 0])
size = 0.3
nodes = mesh.generateMesh([size, size, size], [], [], [], [])

for node in nodes:
    print("Node: ", node.index)
    print("x: ", node.x, " y: ", node.y, " z: ", node.z)
    print("Adj: ", list(node.adj))
