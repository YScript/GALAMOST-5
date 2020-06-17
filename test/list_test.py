# 导入 gamst 包
import gamst

mst=gamst.snapshot.read("example.mst")

# nl=gamst.plists.clist.clist(info=mst, rcut=1.5)
nl=gamst.plists.nlist.nlist(info=mst, rcut=1.5, rbuff=0.1)
nl.calculate()
