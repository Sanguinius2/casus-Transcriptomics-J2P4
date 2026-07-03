# RNA-seq analyse van synoviumweefsel bij reumatoïde artritis


# Inleiding

Reumatoïde artritis (RA) is een chronische auto-immuunziekte waarbij het immuunsysteem lichaamseigen gewrichten aanvalt.
Hierdoor ontstaat ontsteking van het synovium (gewrichtsslijmvlies), wat uiteindelijk kan leiden tot gewrichtsschade.
Hoewel de exacte oorzaak van RA nog niet volledig bekend is, spelen genetische factoren, omgevingsfactoren en ontregeling van het immuunsysteem een belangrijke rol[(Gabriel, 2001; ](https://doi.org/10.1016/S0889-857X(05)70201-5)[Platzer et al., 2019)](https://doi.org/10.1371/journal.pone.0219709).

RA is bekend door complexe veranderingen in genexpressie in het synoviale weefsel, waarbij vooral immuun-gerelateerde pathways sterk geactiveerd zijn. RNA-seq studies hebben aangetoond dat deze veranderingen op transcriptieniveau inzichten kunnen geven in ziekteactiviteit en onderliggende biologische mechanismen[(Platzer et al., 2019)](https://doi.org/10.1371/journal.pone.0219709).

In deze analyse is gebruikgemaakt van RNA-seq data van synoviumbiopten van patiënten met RA en gezonde controles. Het doel was om differentieel tot expressie komende genen en betrokken biologische pathways te identificeren.

Het doel van deze casus is om met behulp van RNA-seq analyse verschillen in genexpressie tussen synoviumweefsel van patiënten met reumatoïde artritis (RA) en gezonde controles te identificeren. Daarnaast wordt onderzocht welke biologische processen en pathways betrokken zijn bij de ziekte door middel van GO enrichment analyse en KEGG pathway analyse.

Voor deze casus zijn de volgende deelvragen opgezet
- 	Welke genen komen differentieel tot expressie in synoviumweefsel van RA-patiënten vergeleken met gezonde controles?
- 	Welke biologische processen zijn oververtegenwoordigd in de differentieel tot expressie komende genen?
- 	Welke immuun-gerelateerde pathways spelen mogelijk een rol bij reumatoïde artritis?

---

## Methoden
De volledige analyse is uitgevoerd met het script [Main_R_File_transcriptomics.R](R_script/Main_R_File_transcriptomics.R). Dit script voert de RNA-seq analyse uit vanaf de alignering van de ruwe FASTQ-bestanden tot en met de differentiële genexpressieanalyse, Gene Ontology (GO) enrichment analyse en KEGG pathway analyse.
De stappen van deze casus bestaan uit de volgende onderdelen(zie ook Figuur 1):
1.	Installeren en laden van de benodigde R-packages.
2.	Opbouwen van een index van het humane referentiegenoom (GRCh38.p14).
3.	Aligneren van paired-end FASTQ-bestanden tegen het referentiegenoom.
4.	Kwantificeren van reads per gen met featureCounts.
5.	Aanmaken van een count matrix.
6.	Uitvoeren van de differentiële genexpressieanalyse met DESeq2.
7.	Visualiseren van de resultaten met een volcano plot.
8.	Uitvoeren van een Gene Ontology enrichment analyse.
9.	Visualiseren van relevante KEGG pathways met pathview.
<p align="center">
  <img src="Figuren/Flowchart.png" alt="Flow" width="600"/>
  
  <em>Figuur 1: Flowchart van uitgevoerde stappen binnen deze casus</em>
</p>


Alle analyses zijn uitgevoerd in R versie 4.6.0 binnen [RStudio](https://docs.posit.co/ide/user/#rstudio-ide-oss-downloads).
De volgende packages zijn gebruikt tijdens de analyse:

*Tabel 1: Alle Packages die in deze casus zijn gebruikt met het versienummer en doel van de packages.*
| Package |	Versie | Doel |
| :------- | :-----: | :----: |
| [Rsubread](https://www.bioconductor.org/packages//release/bioc/html/Rsubread.html) |2.24.0 | Alignering en featureCounts |
| [Rsamtools](https://bioconductor.org/packages/release/bioc/html/Rsamtools.html) | 2.28.0	| Verwerken van BAM-bestanden |
| [DESeq2](https://bioconductor.org/packages/release/bioc/html/DESeq2.html) | 1.50.2	|Differential expression analyse|
| [EnhancedVolcano](https://bioconductor.org/packages/release/bioc/html/EnhancedVolcano.html) | 1.28.2	| Volcano plot |
| [goseq](https://bioconductor.org/packages/release/bioc/html/goseq.html) | 1.62.0	|   GO enrichment analyse |
| [pathview](https://bioconductor.org/packages/release/bioc/html/pathview.html) | 1.70.0	| KEGG pathway visualisatie |
| [KEGGREST](https://bioconductor.org/packages//release/bioc/html/KEGGREST.html) | 1.52.2	| Ophalen KEGG informatie |
| [biomaRt](https://bioconductor.org/packages/release/bioc/html/biomaRt.html) | 2.68.0	| Genannotaties |
| [org.Hs.eg.db](https://bioconductor.org/packages/release/data/annotation/html/org.Hs.eg.db.html) | 3.23.1	| Entrez-ID annotaties |
| [GO.db](https://bioconductor.org/packages/release/data/annotation/html/GO.db.html) | 3.23.1	| Gene Ontology database |
| [tidyverse](https://cran.r-project.org/web/packages/tidyverse/index.html) | 2.0.0 | Dataverwerking en visualisatie |
| [dplyr](https://cran.r-project.org/web/packages/dplyr/index.html) | 1.2.1 | Dataverwerking en visualisatie |
| [ggplot2](https://cran.r-project.org/web/packages/ggplot2/ggplot2.pdf) | 4.0.3 | Dataverwerking en visualisatie |

________________________________________

De analyse is uitgevoerd met R 4.6.0 in het script Main_R_File_transcriptomics.R.
RNA-seq reads van vier controles en vier patiënten met reumatoïde artritis werden uitgelijnd tegen het humane referentiegenoom GCF_000001405.40_GRCh38.p14 van [NCBI](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/) met Rsubread.
Hiervoor werd eerst een genome-index opgebouwd met buildindex(), waarna de FASTQ-bestanden werden gealigneerd met align().
De resulterende .BAM-bestanden werden vervolgens gekwantificeerd met featureCounts, waarbij een count matrix op genniveau werd gegenereerd.
Differentiële genexpressie werd bepaald met DESeq2, waarbij genen als significant werden beschouwd bij een adjusted p-value < 0,05 en een |log2FoldChange| > 1.
De resultaten worden opgeslagen als een CSV-bestand en gevisualiseerd met een volcano plot, gemaakt met EnhancedVolcano.
Vervolgens werd met goseq een Gene Ontology enrichment analyse uitgevoerd.
Hierbij werd gecorrigeerd voor genlengtebias met behulp van biomaRt, waarna de tien meest verrijkte GO Biological Process-termen werden gevisualiseerd.
Ten slotte werden de differentieel tot expressie komende genen geprojecteerd op humane KEGG-pathways met pathview.
Hiervoor werden gen-symbolen omgezet naar Entrez-ID's met org.Hs.eg.db.
De pathwayfiguren en overige resultaten worden automatisch opgeslagen in de repository.


---

## Resultaten

De RNA-seq analyse werd uitgevoerd om verschillen in genexpressie tussen synoviumweefsel van patiënten met reumatoïde artritis (RA) en gezonde controles te identificeren. De analyse bestond uit drie opeenvolgende onderdelen: een differentiële genexpressieanalyse, een Gene Ontology (GO) enrichment analyse en een KEGG pathway analyse.
________________________________________
#Differential expression analyse

Het doel van de differentiële genexpressieanalyse was het identificeren van genen waarvan de expressie significant verschilde tussen RA-patiënten en gezonde controles.
De analyse (Figuur 2) liet duidelijke verschillen in genexpressie zien tussen beide groepen. Er waren 2085 opgereguleerde genen en 2487 neerwaarts gereguleerde genen geïdentificeerd, wat aangeeft dat meerdere biologische processen verschillen tussen gezond weefsel en weefsel afkomstig van patiënten met reumatoïde artritis.


<p align="center">
  <img src="Figuren/VolcanoplotV1.png" alt="Flow" width="600"/>
  
  <em>Figuur 2: Volcano plot van de differentiële genexpressie tussen synoviumweefsel van vier patiënten met reumatoïde artritis en vier gezonde controles. Op de x-as staat de log2FoldChange weergegeven en op de y-as de negatieve log10 van de adjusted p-value. Genen werden als significant beschouwd bij een adjusted p-value kleiner dan 0,05 en een absolute log2FoldChange groter dan 1.(padj)</em>
</p>

________________________________________
#Gene Ontology enrichment analyse

Na het identificeren van de differentieel tot expressie komende genen werd onderzocht welke biologische processen significant oververtegenwoordigd waren. Hiervoor werd een Gene Ontology (GO) enrichment analyse uitgevoerd.
De sterkst verrijkte GO-term was het Immunoglobulin mediated immune response(Figuur 3).
De GO-analyse liet zien dat meerdere immuungerelateerde processen oververtegenwoordigd waren binnen de differentieel tot expressie komende genen.

<p align="center">
  <img src="Figuren/GO-ANALYSE-PLOT-EXTENDED-EDITION" alt="Flow" width="600"/>
  
  <em>Figuur 3: Overzicht van de tien meest significant verrijkte GO Biological Process-termen. Op de x-as is het percentage differentieel tot expressie komende genen binnen iedere GO-term weergegeven. De grootte van de punten geeft het aantal significant veranderde genen binnen iedere categorie weer, terwijl de kleur overeenkomt met de oververtegenwoordigings-p-waarde berekend met goseq.</em>
</p>

________________________________________
#KEGG pathway analyse

Om de biologische betekenis van de differentieel tot expressie komende genen verder te interpreteren, werd een KEGG pathway analyse uitgevoerd.
De pathway analyse liet verhoogde expressie zien van meerdere genen die betrokken zijn bij ontstekingsreacties en adaptieve immuunresponsen. Zowel de pathway Rheumatoid arthritis als de B cell receptor signaling pathway bevatten meerdere genen met verhoogde expressie in de RA-groep ten opzichte van de controles.

<p align="center">
  <img src="Figuren/hsa05323.pathview.png" alt="Flow" width="600"/>
  
  <em>Figuur 4: Visualisatie van de KEGG pathway Rheumatoid arthritis (hsa05323), gegenereerd met pathview. De kleuren geven de richting en grootte van de log2FoldChange weer voor de genen die onderdeel uitmaken van deze pathway.</em>
</p> 
<p align="center">
  <img src="Figuren/hsa04662.pathview.png" alt="Flow" width="600"/>
  
  <em>Figuur 5: Visualisatie van de KEGG pathway B cell receptor signaling pathway (hsa04662). De figuur laat zien welke genen binnen deze signaalroute differentieel tot expressie komen tussen RA-patiënten en gezonde controles.</em>
</p>

---

## Conclusie

Het doel van deze studie was het identificeren van verschillen in genexpressie tussen synoviumweefsel van patiënten met reumatoïde artritis (RA) en gezonde controles met behulp van een RNA-seq analyse.

De differentiële genexpressieanalyse liet zien dat meerdere genen significant verschillend tot expressie kwamen tussen beide groepen.
De Gene Ontology (GO) enrichment analyse toonde aan dat voornamelijk immuun-gerelateerde biologische processen anders waren tussen de twee groepen.
De meest significante GO-term was Immunoglobulin mediated immune response, wat een verhoogde activiteit van B-cellen en antistofgemedieerde afweermechanismen [(Majithia & Geraci, 2007)](https://doi.org/10.1016/j.amjmed.2007.04.005).
Verder liet de KEGG pathway analyse verhoogde activiteit zien binnen de pathways Rheumatoid arthritis (hsa05323) en B cell receptor signaling pathway (hsa04662).

Deze bevindingen beantwoorden de onderzoeksvragen en laten zien dat zowel specifieke genen als immuun-gerelateerde biologische processen en signaalroutes verschillen tussen RA-patiënten en gezonde controles.
De resultaten sluiten goed aan bij de huidige kennis over de pathogenese van reumatoïde artritis. Het is bekend dat B-cellen, auto-antistoffen en chronische ontstekingsprocessen een centrale rol spelen bij het ontstaan en onderhouden van de ziekte [(Smolen et al., 2016)](https://doi.org/10.1056/NEJMra1507093).
Ook eerdere transcriptomicsstudies hebben laten zien dat genexpressie in synoviumweefsel van RA-patiënten wordt gekenmerkt door een verhoogde activiteit van immuun-gerelateerde genen en pathways, waaronder B-celactivatie en antistofgemedieerde immuunresponsen [(Platzer et al., 2019)](https://doi.org/10.1371/journal.pone.0219709). 
De verrijking van de GO-term Immunoglobulin mediated immune response en de gevonden KEGG pathways zijn daarom in overeenstemming met eerder gepubliceerde resultaten. Daarnaast onderstrepen deze bevindingen het belang van B-cellen als potentiële biomarkers en therapeutische aangrijpingspunten bij RA [(Bugatti et al., 2014)](https://doi.org/10.1155/2014/681678).

Een beperking van deze studie is het relatief kleine aantal geanalyseerde samples, waardoor de statistische kracht beperkt is en kleine verschillen mogelijk niet zijn waargenomen.

________________________________________
##Databeheer

Goed databeheer is belangrijk om onderzoek overzichtelijk, reproduceerbaar en betrouwbaar te maken. Tijdens een RNA-seq analyse worden veel verschillende bestanden aangemaakt, zoals FASTQ-bestanden, .BAM-bestanden, count matrices, scripts en figuren.
Door deze bestanden gestructureerd op te slaan blijven analyses eenvoudig terug te vinden en is het mogelijk om de analyses opnieuw uit te voeren en hiermee verder te werken.

De repository is verdeeld in vier hoofdmappen en de README file:

├── README.md

├── Figuren/

├── R_script/

├── Raw Data/

└── Verwerkte Data/

De map Raw Data bevat de originele FASTQ-bestanden. Alle bestanden die tijdens de analyse worden gegenereerd, zoals .BAM-bestanden, de count matrix en DESeq2-resultaten, worden opgeslagen in verwerkte Data. 
De map Figuren bevat alle automatisch gegenereerde figuren en in R File staan het script Main_R_File_transcriptomics.R en het bijbehorende RData-bestand.

Door deze bestanden een duidelijke en consistente naamgeving te geven, bijvoorbeeld Main_R_File_transcriptomics.R, zijn bestanden eenvoudig terug te vinden.
Voor versiebeheer is gebruikgemaakt van GitHub, omdat hier elke versie van de repository te zien is.
Door regelmatig wijzigingen op te slaan met commits blijft de ontwikkelgeschiedenis van het project opgeslagen en kunnen eerdere versies eventueel worden teruggezet of worden ingezien.

GitHub is gebruikt als centrale omgeving voor het beheren van scripts, documentatie en resultaten.
Dankzij een overzichtelijke mappenstructuur, duidelijke bestandsnamen en het centrale script Main_R_File_transcriptomics.R zijn alle onderdelen van de analyse eenvoudig terug te vinden.
Het script Main_R_File_transcriptomics.R bevat commentaarregels waarin iedere stap van de analyse wordt toegelicht.
De README beschrijft de gebruikte dataset, software, workflow en repositorystructuur, zodat andere onderzoekers de analyse eenvoudig kunnen volgen en gebruiken voor eigen onderzoek.
Een ander voordeel van GitHub is dat het open source is, wat inhoud dat iedereen er bij kan. Dit is handig omdat andere onderzoekers het project makkelijk kunnen nadoen en inzien.

Wat wel belangrijk is, is om te denken aan de veiligheid van persoongegevens en andere gevoelige informatie.
De gebruikte RNA-seq data van deze casus bijvoorbeeld zijn afkomstig uit een eerder gepubliceerd onderzoek en bevatten geen herleidbare persoonsgegevens, dit is belangrijk omdat niemand wil dat zijn of haar persoonsgegevens zomaar op straat komen te liggen.
Hierom is het belangrijk om persoonsgegevens te verwijderen of onherkenbaar te maken voordat de code of resultaten online worden gezet.


---

## Bronnen

Bugatti, S., Vitolo, B., Caporali, R., Montecucco, C., & Manzo, A. (2014). *B cells in rheumatoid arthritis: from pathogenic players to disease biomarkers*. BioMed Research International, 2014, 681678. https://doi.org/10.1155/2014/681678

Gabriel, S. E. (2001). *The epidemiology of rheumatoid arthritis*. Rheumatic Disease Clinics of North America, 27(2), 269–281. https://doi.org/10.1016/S0889-857X(05)70201-5

Majithia, V., & Geraci, S. A. (2007). *Rheumatoid Arthritis: Diagnosis and Management*. The American Journal of Medicine, 120(11), 936–939. https://doi.org/10.1016/j.amjmed.2007.04.005

Platzer, A., Nussbaumer, T., Karonitsch, T., Smolen, J. S., & Aletaha, D. (2019). *Analysis of gene expression in rheumatoid arthritis and related conditions offers insights into sex-bias, gene biotypes and co-expression patterns*. PLOS ONE, 14(7). https://doi.org/10.1371/journal.pone.0219709

Smolen, J. S., Aletaha, D., McInnes, I. B. (2016). *Rheumatoid arthritis*. New England Journal of Medicine, 374(21), 2023–2038. https://doi.org/10.1056/NEJMra1507093

Young, M. D., Wakefield, M. J., Smyth, G. K., & Oshlack, A. (2010). *Gene ontology analysis for RNA-seq: accounting for selection bias*. Genome Biology, 11, R14. https://doi.org/10.1186/gb-2010-11-2-r14


