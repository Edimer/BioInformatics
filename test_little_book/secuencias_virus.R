# Genoma completo dengue virus 1 
# https://www.ncbi.nlm.nih.gov/nuccore/NC_001477.1/

library(seqinr)
datos <- read.fasta("sequence_dengue1.fasta")
length(datos$NC_001477.1) #10735
datos$NC_001477.1
class(datos$NC_001477.1)

# Desde R con seqinr
# 1. Apertura de conexión a la base de datos
# 2. consulta
# 3. obtención de secuencias
# 4. cierre de conexión a la base de datos
choosebank("refseqViruses")
dengue <- query(listname = "Dengue1", query = "AC=NC_001477")
seq_dengue <- getSequence(dengue)
anotaciones_dengue <- getAnnot(dengue)
closebank()
length(seq_dengue[[1]]) # 10735


# Secuencia
seq_dengue[[1]]

# Exportar datos fasta
write.fasta(names = "DEN-1", sequences = seq_dengue, file.out = "den1.fasta")
?write.fasta()

# Composición básica de la secuencia
genoma_dengue <- seq_dengue[[1]]
table(genoma_dengue)

# Porcentaje o fracción de G + C
# Contenido GC: https://es.wikipedia.org/wiki/Contenido_GC
sum(table(genoma_dengue)[c(2, 3)]/length(genoma_dengue)) * 100

# con función de seqinr
GC(genoma_dengue)

# Nucleoditos simples: igual a lo anterior con la función count() de seqinr
count(seq = genoma_dengue, wordsize = 1)

# Nucleotidos dobles
count(seq = genoma_dengue, wordsize = 2)

# Nucleotidos triples
count(seq = genoma_dengue, wordsize = 3)

# Contenido de GC local
# En este caso se obtiene gc para cada 2 mil nucleoditos
GC(genoma_dengue[1:1000])
GC(genoma_dengue[1001:2000])
GC(genoma_dengue[2001:3000])
GC(genoma_dengue[3001:4000])
GC(genoma_dengue[4001:5000])
GC(genoma_dengue[5001:6000])
GC(genoma_dengue[6001:7000])
GC(genoma_dengue[7001:8000])
GC(genoma_dengue[8001:9000])
GC(genoma_dengue[9001:10000])
GC(genoma_dengue[10001:length(genoma_dengue)])

# Loop
# La longitud total del genoma es 10735 nucleotidos
# Para analizar la variación local del contenido de GC cada 1000 nucleotidos
# se obtienen "saltos" cada 1000
saltos <- seq(from = 1, to = length(genoma_dengue), by = 2000)
for (i in 1:length(saltos)) {
  if(i < (length(saltos))){
    print(GC(genoma_dengue[saltos[i]:(saltos[i]+1999)]))
  } else {
    print(GC(genoma_dengue[saltos[i]:length(genoma_dengue)]))
  }
}

# Sobrerepresentación y subrepresentación estadística de dinucleotidos
# Rho
# Zscore
# https://www.rdocumentation.org/packages/seqinr/versions/3.6-1/topics/dinucleotides

# Tamaño 1: A y G (frecuencias relativas individuales)
fx_ind <- prop.table(count(seq = genoma_dengue, wordsize = 1))
fx_ind

f_a <- fx_ind["a"]
f_a
f_g <- fx_ind["g"]
f_g

# Tamaño 2:  Frecuencias relativas dobles
fx_doble <- prop.table(count(genoma_dengue, wordsize = 2))
f_ag <- fx_doble["ag"]
f_ag

# Sobre-Sub representación estadística de nucleotidos dobles para AG
f_ag / (f_a * f_g)

# Función Rho
rho(sequence = genoma_dengue, wordsize = 2)


# Ejemplo
# https://www.uniprot.org/uniprot/Q9CD83
# https://www.uniprot.org/uniprot/A0PQ23
library(seqinr)
seq_lepra <- read.fasta("Q9CD83.fasta.txt")
seq_ulcera <- read.fasta("A0PQ23.fasta.txt")

# Desde R
choosebank("swissprot")
lepra <- query("leprae", "AC=Q9CD83")
lepraeseq <- getSequence(lepra)
#ulcera <- query("ulcerans", "AC=A0PQ23") no funcionó
#ulceransseq <- getSequence(ulcerans)
closebank()

# Comparando dos secuencias
sec_ulcera <- toupper(as.character(seq_ulcera$`tr|A0PQ23|A0PQ23_MYCUA`))
sec_lepra <- as.character(lepraeseq[[1]])
length(sec_ulcera)
length(sec_lepra)

# Donde aparee el punto es porque comparten aminoácido
x11()
dotPlot(sec_lepra, sec_ulcera)

# Alineación de secuencias
# Alineanación global o local
# Matriz de sustitución o puntos
library(Biostrings)
