# 导入 gamst 包
import gamst

mst=gamst.snapshot.read("example.mst")

# nl=gamst.plists.clist.clist(info=mst, rcut=1.5)
gg=gamst.chare.particle_set(info=mst, group="all")
print(gg.member)
