TOTAL = 0
CONTADOR = 0

while CONTADOR < 10:
    try:
        nota = float(input(f"{CONTADOR + 1}° Nota: "))
        TOTAL += nota
        CONTADOR += 1
    except ValueError:
        print("Insira um valor numérico.")

print(f"A média das notas é: {TOTAL / 10}")
