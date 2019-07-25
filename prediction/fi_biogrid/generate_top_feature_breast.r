num=100
cutoff=0.01
net=as.matrix(read.table('biogrid_pair.txt'))


fileout = "top_breast.tsv"
f=file('fi_breast.txt','r') # 10 mins
#f=file('fi_breast.txt','r') # 30 mins
line = readLines(f, n = 1)
id_all = unlist(strsplit(line,split='\t'))[-1]
num0 = length(id_all)

p=NULL

i=0
while (TRUE) {
    print(i)
    i=i+1
    line = readLines(f, n = 1)
    if ( length(line) == 0 ) { break }
    tmp = unlist(strsplit(line,split='\t'))
    the_id = unlist(strsplit(tmp[1],split='\\.'))[1]
    val = as.numeric(tmp[-1])

    list1 = id_all[order(val,decreasing=T)[1:num]]
    list2 = unique(as.vector(net[net[,1]==the_id | net[,2]==the_id,]))
    list2 = intersect(list2, id_all)
    list_overlap = intersect(list1, list2)

    num1 = length(intersect(list1,list2))

    d = data.frame(gene.in.interest=c(num1, num-num1), gene.not.interest=c(length(list2)-num1, num0-num-length(list2)+num1))
    the_p = fisher.test(d)$p.value
    p = c(p,the_p)
    if (the_p < cutoff){
        write(paste0(the_p, '\t', paste0(list_overlap, collapse=','), '\t', paste0(list1, collapse=',')), file=fileout, append=TRUE)
    }
}
close(f)


