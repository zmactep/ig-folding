import sys
import heapq


class Entry:
    energy = 0
    model = ''

    def __init__(self, energy, model):
        self.energy = energy
        self.model = model

    def __lt__(self, other):
        return self.energy < other.energy

    def __ge__(self, other):
        return self.energy >= other.energy


def main():
    file_name = sys.argv[1]
    top_n = int(sys.argv[2])
    heap = [Entry(0, '')] * top_n

    with open(file_name) as f:
        next(f)
        next(f)
        for line in f:
            x = line.split()
            e = Entry(-float(x[4]), x[-1])
            if e >= heap[0]:
                heapq.heapreplace(heap, e)
    heap.sort(reverse=True)
    for i in range(top_n):
        print("Best " + str(i+1) + " model: " + heap[i].model + " with score(I_sc): " + str(-heap[i].energy))


if __name__ == "__main__": main()