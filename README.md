# HAA Analysis – Four Kaon Final State

This repository contains the analysis framework for a search for light pseudoscalars (a) produced in Higgs boson decays (H → aa → 4K), specifically targeting the four-kaon final state.

## 📁 Project Overview

- **Goal**: Search for exotic Higgs decays involving pseudoscalars using CMS data
- **Signal**: H → aa → (K⁺K⁻)(K⁺K⁻)
- **Data**: CMS Run 2 Data / private samples
- **Tools**: ROOT, Combine, CMSSW, KingMaker

## 🛠️ Setup Instructions

### 1. Clone KingMaker

```bash
git clone git@github.com:KIT-CMS/KingMaker.git
cd KingMaker
source setup.sh
cd ..
```

### 2. Clone the repository

```bash
git clone git@github.com:jhornung97/haa_analysis.git
cd haa_analysis
cp -r haa_analysis/KingMaker/* KingMaker 
```
