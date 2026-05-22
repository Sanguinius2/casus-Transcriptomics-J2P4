# casus
<p align="center">
  <img src="Figuren/hsa04662.png" alt="Flow" width="600"/>
</p>


# Casus Transcriptomis

Welkom in **Wizarding Spell Metrics**, waar de data nep is, maar de structuur *verrassend oké* is. Deze repo is gemaakt om te laten zie hoe je een bio-informatica project structureert in GitHub met behulp van een onzin dataset. files mogen gedownload en hergebruikt worden (zoals [deze README.md file](README.md)). Vraag ook gerust AI voor tips!

Perfect voor:
- Wegwijs worden in GitHub
- Het leren structureren van data

> `You Shall Not PASS`

## 📁 Inhoud mapjes structuur

- `data/raw/` –  
- `data/processed` - verwerkte datasets gegenereerd met scripts 
- `R Files/` – R Files om de 
- `resultaten/` - grafieken en tabellen
- `bronnen/` - gebruikte bronnen 
- `README.md` - het document om de tekst hier te genereren
- `data_Beheer/` - de info over github en waarom het belangrijk is 

---

## Introductie

Deze GitHub pagina gaat over Rheuma. Met behulp van de R files kan je zelf ook *UITVOEREN!!!!*.

Als je naar het menselijk genoom wil, [klik hier](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/) om naar .

## Methoden

De data is verstrekt door de Afdeling Magische Wetshandhaving en Ollivanders Wandwinkel Archieven. 

De ruwe data van spreuken is eerst bewerkt voor analyse met behulp van [scripts/01_clean_spell_data.R](scripts/01_clean_spell_data.R). Vervolgens zijn de spreuken geanalyseerd op kracht en nauwkeurigheid met [scripts/02_spell_analysis.R](scripts/02_spell_analysis.R).

## 📊 Resultaten

Om inzicht te krijgen in eigenschappen van de te gebruiken spreuken is er een overzicht gemaakt, te vinden in [deze tabel](resultaten/top_10_spells.csv). Onvergeeflijke vloeken zijn niet meegenomen in dit overzicht. 

Om een afweging te maken welke spreuken het meest effectief zijn, is er onderzocht of er een verband te vinden is tussen kracht en accuraatheid. In [het resultaat hiervan](resultaten/spell_power_vs_accuracy.png) is te zien dat er een negatieve daling lijkt te zijn in kracht als de accuraatheid toeneemt. Een uitschieter is de onvergeeflijke vloek *Avada Kedavra*, met zowel hoge kracht als accuraatheid. 

## Conclusie

Spreuken met meer accuraatheid lijken minder krachtig te zijn. Een uitzondering op deze trend is de onvergeeflijke vloek *Avada Kedavra*, welke beter niet gebruikt kan worden. 





