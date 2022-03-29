quantFiles=(`ls salmon/mRNA/*/quant.sf`)
 
 mkdir -p quantFiles

 for f in ${quantFiles[@]}
 do
 	newName=${f#salmon/mRNA/}
 	newName=${newName%/quant.sf}
 	cp $f ./quantFiles/${newName}_quant.sf
 done

