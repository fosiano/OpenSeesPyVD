#annullo alcune variabili per evitare che contengono dati di esecuzioni non andate a buon fine 
#(utile in fase di test)
in_file=None
out_file=None
line=None
dati = None
#fine

#memorizzazione del nome dei file di input/output a meno dell'estensione
nomefile="AccTH-MSS-SLV"
for i in range(1,11):
    nfile=nomefile+" ("+str(i)+")"
    #fine
    
    #apertura file di input in modalità di lettura
    fileIN=nfile+".txt" #nome file di input
    in_file = open(fileIN,"r")
    #fine
    
    #apertura file di output in modalità di scrittura
    fileOUT=nfile+".th" #nome file di output
    out_file = open(fileOUT,"w")
    #fine
    
    #echo
    print ("Lettura file in corso " + fileIN + "....")
    print ("Scrittura file in corso " + fileOUT + "....")
    print ("Echo ... ")
    #echo
    
    #lettura due righe a vuoto per presenza intestazioni
    line=in_file.readline()
    line=in_file.readline()
    riga=2
    #fine
    
    #lettura prima riga utile
    line=in_file.readline()
    #fine
    
    #ciclo fino alla fine del file di input
    amax=0.0
    while (not (line == "EOF\n" or line=="")):
        #memorizzazione della riga corrente incrementando di una unità
        riga+=1
        #fine
        
        #interpretazione dell'intera riga di testo in valori elementari
        dati = line.split() #i valori sono separati da " "
        #fine
        col2_string = dati[1]
        acc=float(col2_string)/9.81
        if abs(acc) > amax:
            amax=abs(acc)
        out_file.write(str(acc)+"\n") 
        line=in_file.readline() 
    #fine ciclo per raggiungimento della fine del file di input
    
    #chiusura file aperti in lettura/scrittura
    in_file.close()
    out_file.close()
    #fine
    
    #echo
    print("\n")
    print ("Lettura file " + fileIN + " COMPLETATA")
    print ("Scrittura file " + fileOUT + " COMPLETATA")
    print(amax)
    #echo





