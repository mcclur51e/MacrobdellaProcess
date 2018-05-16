library(RAM)

OTU.MdILF = taxdt.MacILF$Number
OTU.MdInt = taxdt.MacInt$Number
OTU.MdBlad = taxdt.MacBlad$Number


x = list(
  "ILF" = OTU.MdILF,
  "Bladder" = OTU.MdBlad,
  "Intestinum" = OTU.MdInt
)

y = list(
  "ILF" = OTU.MdILF,
  "Bladder" = OTU.MdBlad
)

group.venn(vectors=x, label=FALSE)
ggsave(grid.draw(RAM::group.venn(vectors=x, label=FALSE)), filename="vennCore.png", width=12,height=8)

group.venn(vectors=y)
ggsave(grid.draw(RAM::group.venn(vectors=y)), filename="vennCore_listBladILF.png", width=12,height=8)

#library(VennDiagram)
#ggsave(grid.draw(VennDiagram::venn.diagram(x, NULL)), filename="vennCore.png", width=12,height=8)



