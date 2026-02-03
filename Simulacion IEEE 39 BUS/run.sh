# Este script automatiza la ejecucion del modelo main.gms para la generacion de archivos csv
	
# Se especifican las entradas y directorios de salida
	
# El script cumple la misma funcion que run.ps1
# Se ejecuta en git bash (linux)
	
# Para ejecutar este script, escribir en la terminal: ./run.sh

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
