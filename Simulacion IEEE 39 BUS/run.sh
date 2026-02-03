for i in {0..47}
do
  echo "Caso $i"
  gams main.gms \
    --DATA=Inputs/Input_${i}.DAT \
    --PRICE=Precios/Precios_${i}.csv \
    --HIST=Historia/Historia_${i}.csv \
    --GEN=Generacion/Generacion_${i}.csv\
    lo=2
done
