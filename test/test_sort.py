import gamst
	
mst = gamst.snapshot.read("lj.mst")

fn = gamst.tinker.sort(info=mst, period=1)
fn.data.calculate(0)

