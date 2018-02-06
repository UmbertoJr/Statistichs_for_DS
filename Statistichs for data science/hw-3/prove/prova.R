install.packages("igraph")
library(igraph)
library(Matrix)

g <- graph.empty(n=10,directed = T)
plot(g)
h <- graph.full(n=5,directed = F,loops = T)
plot(h)

l <- graph.star(n=10, mode='mutual')
plot(l)

h <- graph.ring(n=10)
plot(h)

edges = c(1,2,3,2,2,4)
m <- graph(edges, n=max(edges))
plot(m)

vcount(m)
ecount(m)
neighbors(m, V(m)[2], mode=3)

incident(m,V(m)[2], mode = c("all", "out", "in", "total"))
incident(m,V(m)[2], mode ="total")

get.edgelist(m)

is.directed(m)
are.connected(m, 1, 4)



mat <-get.adjacency(h)
plot(m)

V(m)
E(m,P=NULL)

?add.vertices
?add.edges

plot(m)
V(m)

gi <- make_empty_graph() %>%
  add_vertices(3, color = "red") %>%
  add_vertices(2, color = "green") %>%
  add_edges(c(1,2, 2,3, 3,4, 4,5))
gi
V(gi)[[]]
plot(gi)


gr <- make_empty_graph(n = 5) %>%
  add_edges(c(1,2, 2,3, 3,4, 4,5)) %>%
  set_edge_attr("color", value = "red") %>%
  add_edges(c(5,1), color = "green")
E(gr)[[]]
plot(gr)

--------------------

  
n <- sample(5:10,1)
gra <- graph.ring(n)
vcount(gra)
ecount(gra)

neighbors(gra,5)
incident(gra,5)
are.connected(gra,1,3)
are.connected(gra,3,4)

plot(gra, layout= layout.fruchterman.reingold.grid, vertex.label=V(gra)$number, edge.arrow.size=0.5)

----------------------------

ph <- graph.empty(5, directed = F)    
new <- c(1,3,1,5,2,5,4,5)
ph <- add.edges(ph,new)
plot(ph)

ph <- add.vertices(ph,1)
ph <- add.edges(ph,c(6,5))
V(ph)$name <- letters[1:vcount(ph)]

# max(degree(ph))
E(ph)$pesi <- runif(ecount(ph))
get.adjacency(ph, attr = 'pesi')
get.adjlist(ph)

plot(ph)
write.graph(ph,'grafo1.dl','pajek')
write.graph(ph, 'grafo1_2.txt', 'edgelist')

advice_data_frame <- read.table('http://sna.stanford.edu/sna_R_labs/data/Krack-High-Tec-edgelist-Advice.txt')
big_g<- graph.data.frame(advice_data_frame)
length(E(big_g))
V(big_g)
?graph.data.frame

la <- graph.star(8)
num <- sample(1:50,vcount(la), replace = T)
V(la)$number <- num
V(la)$col <- 'grey'
V(la)[number<30]$col <- 'green'
plot(la,layout=layout.circle, vertex.color=V(la)$col, vertex.label=V(la)$number, edge.arrow.size=0.09)


bo <- graph.full(10)
E(bo)$pes <- runif(ecount(bo))
E(bo)$lar = 0.5
E(bo)$col = 'red'
E(bo)[pes<0.5]$lar=1.5
E(bo)[pes<0.5]$col = 'green'
plot(bo, layout=layout.circle,edge.width=E(bo)$lar, edge.color=E(bo)$col)
E(bo)$pes
