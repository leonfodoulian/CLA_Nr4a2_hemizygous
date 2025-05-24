# Cluster attributes
cluster.attributes <- list(
  # All cells without filtering
  "all_cells" = list(
    broad.clusters = list(
      # Cluster names
      cluster.names = c(
        "1" = "microglia",
        "2" = "microglia",
        "3" = "vascular leptomeningeal cell",
        "4" = "layer 6 glutamatergic neuron",
        "5" = "astrocyte",
        "6" = "shell-like glutamatergic neuron",
        "7" = "striatal medium spiny neuron",
        "8" = "endothelial cell",
        "9" = "microglia",
        "10" = "oligodendrocyte progenitor cell",
        "11" = "noise",
        "12" = "claustrum projection neuron",
        "13" = "oligodendrocyte",
        "14" = "interneuron",
        "15" = "fibroblast"
      ),
      
      # Cluster abbreviations
      cluster.abbreviations = c(
        "1" = "micro",
        "2" = "micro",
        "3" = "VLMC",
        "4" = "L6",
        "5" = "astro",
        "6" = "shell-like",
        "7" = "MSN",
        "8" = "endo",
        "9" = "micro",
        "10" = "OPC",
        "11" = "noise",
        "12" = "CLA",
        "13" = "oligo",
        "14" = "IN",
        "15" = "fibro"
      ),
      
      # Cluster breaks
      cluster.breaks = c("CLA", "shell-like", "L6", "IN", "MSN", "astro", "OPC", "oligo", "micro", "endo", "VLMC", "fibro", "noise"),
      
      # Cluster colors
      cluster.colors = c(
        "CLA" = "#4B0082",
        "shell-like" = "#F1790F",
        "L6" = "#42AED4",
        "IN" = "#278081",
        "MSN" = "#91C548",
        "astro" = "#CBC597",
        "OPC" = "#A98253",
        "oligo" = "#7D8D87",
        "micro" = "#92876C",
        "endo" = "#79735B",
        "VLMC" = "#655E4B",
        "fibro" = "#504A3D",
        "noise" = "#D3D3D3"
      )
    ),
    
    # Clustering result
    clustering.result = "SCT_snn_res.0.3"
  ),
  
  # Neurons without filtering
  "neurons_unfiltered" = list(
    broad.clusters = list(
      # Cluster names
      cluster.names = c(
        "1" = "layer 6a glutamatergic neuron",
        "2" = "shell projection neuron",
        "3" = "claustrum projection neuron",
        "4" = "striatal indirect pathway medium spiny neuron",
        "5" = "striatal direct pathway medium spiny neuron",
        "6" = "layer 2/3, layer 5b, piriform glutamatergic neuron",
        "7" = "noise",
        "8" = "layer 6b glutamatergic neuron",
        "9" =  "interneuron",
        "10" = "striosome medium spiny neuron",
        "11" = "microglia"
      ),
      
      # Cluster abbreviations
      cluster.abbreviations = c(
        "1" = "L6a",
        "2" = "shell",
        "3" = "CLA",
        "4" = "iMSN",
        "5" = "dMSN",
        "6" = "L2/3-L5b-pir",
        "7" = "noise",
        "8" = "L6b",
        "9" = "IN",
        "10" = "sMSN",
        "11" = "micro"
      ),
      
      # Cluster breaks
      cluster.breaks = c("CLA", "shell", "L6a", "L6b", "L2/3-L5b-pir", "IN", "dMSN", "iMSN", "sMSN", "micro", "noise"),
      
      # Cluster colors
      cluster.colors = c(
        "CLA" = "#4B0082",
        "shell" = "#FF4D00",
        "L6a" = "#1E90FF",
        "L6b" = "#66CDAA",
        "L2/3-L5b-pir" = "#E4A51F",
        "IN" = "#278081",
        "dMSN" = "#7DAE39",
        "iMSN" = "#6B982D",
        "sMSN" = "#50761D",
        "micro" = "#92876C",
        "noise" = "#D3D3D3"
      )
    ),
    
    # Clustering result
    clustering.result = "SCT_snn_res.0.2"
  ),
  
  # Neurons filtered
  "neurons_filtered" = list(
    broad.clusters = list(
      # Cluster names
      cluster.names = c(
        "1" = "layer 6a glutamatergic neuron",
        "2" = "claustrum projection neuron",
        "3" = "striatal indirect pathway medium spiny neuron",
        "4" = "shell projection neuron",
        "5" = "striatal direct pathway medium spiny neuron",
        "6" = "layer 6a glutamatergic neuron",
        "7" = "shell projection neuron",
        "8" = "layer 6a glutamatergic neuron",
        "9" = "layer 6b glutamatergic neuron",
        "10" = "layer 2/3 glutamatergic neuron",
        "11" = "interneuron",
        "12" = "layer 6a glutamatergic neuron",
        "13" = "striosome medium spiny neuron",
        "14" = "piriform glutamatergic neuron",
        "15" = "shell projection neuron",
        "16" = "mixed striatal medium spiny neuron",
        "17" = "claustrum projection neuron",
        "18" = "layer 5b glutamatergic neuron"
      ),
      
      # Cluster abbreviations
      cluster.abbreviations = c(
        "1" = "L6a",
        "2" = "CLA",
        "3" = "iMSN",
        "4" = "shell",
        "5" = "dMSN",
        "6" = "L6a",
        "7" = "shell",
        "8" = "L6a",
        "9" = "L6b",
        "10" = "L2/3",
        "11" = "IN",
        "12" = "L6a",
        "13" = "sMSN",
        "14" = "pir",
        "15" = "shell",
        "16" = "mMSN",
        "17" = "CLA",
        "18" = "L5b"
      ),
      
      # Cluster breaks
      cluster.breaks = c("CLA", "shell", "L6a", "L6b", "L2/3", "L5b", "pir",
                         "IN", "dMSN", "iMSN", "mMSN", "sMSN"),
      
      # Cluster colors
      cluster.colors = c(
        "CLA" = "#4B0082",
        "shell" = "#FF4D00",
        "L6a" = "#1E90FF",
        "L6b" = "#66CDAA",
        "L2/3" = "#E2CC3B",
        "L5b" = "#F0880E",
        "pir" = "#E0A01A",
        "IN" = "#278081",
        "dMSN" = "#7DAE39",
        "iMSN" = "#6B982D",
        "mMSN" = "#A9D46A",
        "sMSN" = "#50761D"
      )
    ),
    
    subclusters = list(
      # Cluster names
      cluster.names = c(
        "1" = "layer 6a glutamatergic neuron 1",
        "2" = "claustrum projection neuron 1",
        "3" = "striatal indirect pathway medium spiny neuron",
        "4" = "shell projection neuron 1",
        "5" = "striatal direct pathway medium spiny neuron",
        "6" = "layer 6a glutamatergic neuron 2",
        "7" = "shell projection neuron 2",
        "8" = "layer 6a glutamatergic neuron 3",
        "9" = "layer 6b glutamatergic neuron",
        "10" = "layer 2/3 glutamatergic neuron",
        "11" = "interneuron",
        "12" = "layer 6a glutamatergic neuron 4",
        "13" = "striosome medium spiny neuron",
        "14" = "piriform glutamatergic neuron",
        "15" = "shell projection neuron 3",
        "16" = "mixed striatal medium spiny neuron",
        "17" = "claustrum projection neuron 2",
        "18" = "layer 5b glutamatergic neuron"
      ),
      
      # Cluster abbreviations
      cluster.abbreviations = c(
        "1" = "L6a 1",
        "2" = "CLA 1",
        "3" = "iMSN",
        "4" = "shell 1",
        "5" = "dMSN",
        "6" = "L6a 2",
        "7" = "shell 2",
        "8" = "L6a 3",
        "9" = "L6b",
        "10" = "L2/3",
        "11" = "IN",
        "12" = "L6a 4",
        "13" = "sMSN",
        "14" = "pir",
        "15" = "shell 3",
        "16" = "mMSN",
        "17" = "CLA 2",
        "18" = "L5b"
      ),
      
      # Cluster breaks
      cluster.breaks = c("CLA 1", "CLA 2", "shell 1", "shell 2", "shell 3",
                         "L6a 1", "L6a 2", "L6a 3", "L6a 4", "L6b", "L2/3", "L5b", "pir",
                         "IN", "dMSN", "iMSN", "mMSN", "sMSN"),
      
      # Cluster colors
      cluster.colors = c(
        "CLA 1" = "#8A4FBA",
        "CLA 2" = "#2D0050",
        "shell 1" = "#FF9466",
        "shell 2" = "#FF4D00",
        "shell 3" = "#A02900",
        "L6a 1" = "#7DBEFE",
        "L6a 2" = "#55A0EA",
        "L6a 3" = "#3E74AD",
        "L6a 4" = "#155498",
        "L6b" = "#66CDAA",
        "L2/3" = "#E2CC3B",
        "L5b" = "#F0880E",
        "pir" = "#E0A01A",
        "IN" = "#278081",
        "dMSN" = "#7DAE39",
        "iMSN" = "#6B982D",
        "mMSN" = "#A9D46A",
        "sMSN" = "#50761D"
      )
    ),
    # Clustering result
    clustering.result = "SCT_snn_res.1.2"
  )
)
