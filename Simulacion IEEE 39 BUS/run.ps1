<#
- Este script automatiza la ejecucion del modelo main.gms para la generacion de archivos csv
	
- Se especifican las entradas y directorios de salida
	
- El script cumple la misma funcion que run.sh
	A diferencia de la version .sh, esta se ejecuta desde powershell (nativo de windows)
	
- Para ejecutar este script, escribir en la terminal: run.ps1 o .\run.ps1

- En caso aun no se puedan ejecutar scripts .ps1, escribir:
Set-ExecutionPolicy -Scope CurrentUser -ExecutionPolicy RemoteSigned
#>

0..47 | ForEach-Object {
    Write-Host "Caso $_"

    gams main.gms `
        --DATA=Inputs/Input_$_.DAT `
        --PRICE=Precios/Precios_$_.csv `
        --HIST=Historia/Historia_$_.csv `
        --GEN=Generacion/Generacion_$_.csv `
        lo=2
}

