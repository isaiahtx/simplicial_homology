import json
from smithnormalform import matrix, snfproblem, z
import time

class up_tree:
    def __init__(self,size):
        self.size = size
        self.nodes = [-1] * size

    def find(self,x):
        if self.nodes[x] < 0: return x
        if self.nodes[self.nodes[x]] < 0: return self.nodes[x]
        self.nodes[x] = self.find(self.nodes[x])
        return self.nodes[x]

    def union(self,x,y):
        sentinel_x = self.find(x)
        sentinel_y = self.find(y)
        if sentinel_x == sentinel_y: return
        if self.nodes[sentinel_y] <= self.nodes[sentinel_x]:
            self.nodes[sentinel_y] += self.nodes[sentinel_x]
            self.nodes[sentinel_x] = sentinel_y
        else:
            self.nodes[sentinel_x] += self.nodes[sentinel_y]
            self.nodes[sentinel_y] = sentinel_x
    
    def components(self):
        roots = set({})
        for i in range(self.size):
            roots.add(self.find(i))
        return len(roots)

def big_to_small(fset):
    while True:
        additions_made = False
        for i in range(len(fset)):
            for j in range(len(fset[i])):
                f = fset[i][:j] + fset[i][j+1:]
                if f not in fset and len(f) > 0:
                    fset.append(f)
                    additions_made = True
        if not additions_made: break

    return fset

def main():
    print("Line segment disjoint union (S^1 v S^2):")
    simp_hom(big_to_small([[0,1,2],[0,1,3],[0,2,3],[1,2,3],[3,4],[1,4],[5,6]]))
    print()

    print("RP^2:")
    simp_hom(big_to_small([[0,4,5], [0,1,5], [1,3,5], [0,3,4], [0,2,3], [0,1,2], [1,2,4], [1,3,4], [2,4,5], [2,3,5]]))
    print()

    print("S^1:")
    simp_hom(big_to_small([[0,1],[1,2],[2,0]]))
    print()
    
    print("Klein bottle:")
    simp_hom(big_to_small([[0,1,6],[1,6,7],[1,2,7],[2,7,8],[0,2,8],[0,3,8],[3,6,7],[3,4,7],[4,7,8],[4,5,8],[3,5,8],[3,5,6],[0,3,4],[0,1,4],[1,4,5],[1,2,5],[2,5,6],[0,2,6]]))


def get_nontrivial_elementary_divisors(mat):
    if mat == [[]]: return [],0
    rows = len(mat)
    cols = len(mat[0])
    smat = []
    for i in range(rows):
        for j in range(cols):
            smat.append(z.Z(mat[i][j]))
    smat = matrix.Matrix(rows,cols,smat)
    snfmat = snfproblem.SNFProblem(smat)
    while True:
        try:
            snfmat.computeSNF()
            break
        except:
            time.sleep(0.2)
            continue
    elem_divisiors = []
    r = 0
    for i in range(min(rows,cols)):
        divisor = int(str(snfmat.J.get(i,i)))
        if divisor != 0:
            r += 1
        if not (-2 < divisor < 2):
            elem_divisiors.append(int(str(snfmat.J.get(i,i))))
    return elem_divisiors, r

def get_matrix(s1,s2,o1,o2):
    rank_s1 = len(s1)
    rank_s2 = len(s2)
    if rank_s1 == 0:
        return [[]]
    mat = []
    for i in range(rank_s2):
        mat.append([0] * rank_s1)
    for i in range(len(s1)):
        boundary = differential(s1[i])
        for j in range(len(boundary)):
            mat[(o2[json.dumps(boundary[j][1])])][o1[json.dumps(s1[i])]] = boundary[j][0]

    return mat

def differential(gen):
    # calculates only differentials of generators.
    output = []
    c = -1
    for i in range(len(gen)):
        c *= -1
        output.append((c,gen[:i]+gen[i+1:]))
    return output

def simp_hom(fset):
    biggest_simplex = -1
    for i in fset:
        if len(i) > biggest_simplex: biggest_simplex = len(i)

    orders = []
    scomplex = []
    for i in range(biggest_simplex+1):
        scomplex.append([])
        orders.append(dict({}))
    for i in fset:
        if len(i) > 0: scomplex[len(i)-1].append(i)


    for i in range(1,len(scomplex)):
        for j in range(len(scomplex[i])):
            scomplex[i][j].sort()
        scomplex[i].sort()

    for i in range(len(scomplex)):
        for j in range(len(scomplex[i])):
            orders[i][json.dumps(scomplex[i][j])] = j
    
    divisors_and_ranks = [(0,[])]

    for i in range(biggest_simplex):
        mat = get_matrix(scomplex[i+1],scomplex[i],orders[i+1],orders[i])
        ntemd,r = get_nontrivial_elementary_divisors(mat)
        divisors_and_ranks.append((r,ntemd))

    homology_groups = []

    for i in range(biggest_simplex):
        m = len(scomplex[i])
        s = divisors_and_ranks[i][0]
        r = divisors_and_ranks[i+1][0]
        homology_groups.append([m-r-s] + divisors_and_ranks[i+1][1])

    for i in range(len(homology_groups)):
        print('H' + str(i) + ': ',end='')
        if homology_groups[i] == [0]:
            print(0,end='')
        else:
            print("Z",end='')
            if homology_groups[i][0] not in [0,1]:
                print('^' + str(homology_groups[i][0]),end='')
        if len(homology_groups[i]) > 1:
            for j in range(1,len(homology_groups[i])):
                print(' ⊕ Z/' + str(homology_groups[i][j]) + 'Z',end='')
        print()
    print("Hn: 0\t for all n>"+str(i))




if __name__ == "__main__":
    main()
