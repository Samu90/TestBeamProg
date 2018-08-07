# TestBeamProg

-Wiki base delle ntuple utilizzate per l'analisi dei dati del Test Beam di MAGGIO 2018 A T10

-plotWF.C(Conf,nevento) plotta le 6 forme d'onda acquisite di un evento (con 0 le plotta in sequenza, non si ferma)

-ploWF_graph.C(Conf) plotta gli istogrammi delle ampiezze massime e crea un file histoConf.root con gli istogrammi

-plotHisto.C(histoConf1,histoconf2) plotta gli istogrammi delle ampiezze massime delle forme d'onda

-plotWF_fit.C(Conf) plotta gli istogrammi delle ampiezze massime delle forme d'onda e fa un fit con la funzione di Landau

-plotWF_cut.C(Conf) plotta gli istogrammi delle ampiezze massime delle forme d'onda e fa un fit con la funzione di Landau e taglia gli eventi in un intervallo (0.8*MIP,3*MIP)

-plotWF_tamp.C(Conf) plota gli istogrammi 2D di t_r-t_MCP t_l-t_MCP t_ave-t_MCP in funzione dell'ampiezza massima, fa un fit sui valori medi delle variabili temporali per ogni valore di dei bin dell'ampiezza massima.

-plot_lsign.C(Conf) plotta una gaussiana con la sigma dei tempi calcolata senza correggere con il time walk.

-plotWF_su.C(Conf) plotta i valori di amp_max in funzione di t_left-t_right

-plotWF_corr.C(Conf) plotta i grafici di _tamp.C insime a quelli corretti con l'amplitude walk, sono stati usati i tempi calcolati con soglia LED30 (30 u.d.). Plotta infine la gaussiana con la distribuzione 