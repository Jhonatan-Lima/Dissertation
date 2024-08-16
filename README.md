
## Apresentação

Este repositório é um complemento do trabalho de pesquisa intitulado
**“Agregação espacial das populações de peixes da planíce de inundação
do Alto rio Paraná”**, qual realizei durante minha dissertação de
mestrado no Programa de Ecologia e Evolução da Universidade Federal de
Goiás sobre a orientação do **Dr. Luis Mauricio Bini** e coorientação da
**Dra. Rafaela Vendrametto Granzotti**.

### Organização do repositório

Os arquivos estão distribuídos em pastas, a organização das pastas busca
ser o mais intuitivo possível de forma a facilitar a navegação no
repositório. Abaixo é possível encontrar uma breve descrição de cada
pasta e dos principais arquivos de dados, todos os arquivos não
descritos podem ser gerados através dos scripts.

#### Pasta `data`

A pasta `data` contém todos os dados nescessários para realizar as
análises.  

- `cpue_peld.txt` : é uma matriz da cpue com espécies nas colunas e
  coletas nas linhas.  
- `peld_complete.txt` : são os dados tabulares onde cada linha
  representa um indivíduo coletado.  
- `species_traits.csv` : é uma matriz com características das espécies
  nas colunas e espécies nas linhas.  
- `mean_hydro.txt` : é uma matriz com a média mensal do nível
  hidrológico do rio Paraná.  
- `fish_phylo.nex` : é uma filogenia das espécies de peixes (os métodos
  do segundo capítulo detalha sua construção).  

**OBS**  
Os arquivos `cpue_peld.txt` e `peld_complete.txt` não estarão disponível
publicamente neste repositório.  

#### Pasta `doc`

A pasta `doc` contém a versão atualizada do manuscrito.

#### Pasta `output`

Na pasta `output` são encontrados as tabelas e figuras inseridas no
manuscrito.

#### Pasta `R`

A pasta `R` contém todos os **scripts** nescessários para ler dados e
realizar as análises deste trabalho, abaixo explico a lógica utilizada
para nomeação dos arquivos.  

- `Script + 00` : os arquivos com **“script”** no início do nome contém
  o código principal para realizar as análises deste trabalho, e o
  número informando logo após indica a sequência que a análise deve ser
  realizada.  
- `Letra maiúscula` : cada letra informa o tipo de análise que o script
  realiza.  
  - **C** = leitura e limpeza/processamento de dados  
  - **D** = análises  
  - **S** = gráficos/figuras  
- `Descrição` : após cada letra terá uma descrição sobre o que o script
  faz especificamente.  
- `Functions` : é uma pasta com arquivos que contém uma função
  necessária para rodar o **“script principal”**.

### Reprodutibilidade das análises

Este repositório também conta com sistema de manejo de versões de
pacotes por meio do pacote {`renv`}, qual garantirá o funcionamento
correto de todas as funções que são executadas nos scripts. Para isso
basta executar o código abaixo:

**Importante**  
O código fará a instalação automática em seu computador de todos os
pacotes nescessários para rodar as análises do projeto na mesma versão
que foi utilizada durante a elaboração dos scripts.

``` r
# Install package "renv"
install.packages("renv")

# Install versions of the packages used
renv::restore()
```

### Contato

Se você encontrar algum bug ou tiver alguma dúvida, pode entrar em
contato comigo através de <lima.fishing@gmail.com>.
