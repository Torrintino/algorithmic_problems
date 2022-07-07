import matplotlib.pyplot as plot
import numpy


tries = list(numpy.fromfile("tries.bin", dtype=numpy.uint32))
clause_size = [i for i in range(500, 2100, 100)]

for t in tries:
    actual = list(numpy.fromfile("actual{}.bin".format(t), dtype=numpy.uint32))
    expected = list(
        numpy.fromfile("expected{}.bin".format(t), dtype=numpy.uint32))
    fig, ax = plot.subplots()
    ax.plot(clause_size, expected, label="Expected at least")
    ax.plot(clause_size, actual, label="Actual")
    ax.plot(clause_size, clause_size, label="Maximum possible")
    ax.set_xlabel("Number of Clauses in Formula")
    ax.set_ylabel("Satisfied Clauses")
    ax.set_title(
        "Average satisfied clauses in random MAX-SAT with {} assignments".format(t))
    ax.legend()
