# Cálculo de Órbitas

Neste trabalho, usaremos o método numérico de Runge-Kutta de 4ª Ordem para simularmos,
computacionalmente, a órbita de um objeto ao redor de um buraco negro, e ainda a trajetória da
luz que chega a um observador distante, segundo a métrica de Schwarzschild. O código é capaz de calcular trajetórias orbitais e gerar dados para análise, como posição e velocidade em coordenadas polares.

## Funcionalidades

O projeto oferece as seguintes funcionalidades:


1. **Simulação das Órbitas**:
   - `orbit_data(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)`: Simula a órbita usando o método de Runge-Kutta de 4ª ordem. Os dados são salvos em um arquivo chamado `orbit_data.txt`.
   - `ver_orbitas(h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi)`: Simula a órbita da luz usando o método de Runge-Kutta de 4ª ordem. Os dados são salvos em um arquivo chamado`eye_data.txt`.

2. **Geração de Dados**:
   - Os resultados das simulações são salvos em arquivos de texto (`orbit_data.txt` e `eye_data.txt`), contendo informações sobre o tempo (`t`), ângulo (`phi`), e distância radial (`r`).

3. **Plotando os Dados**
   - Os dados gerados podem ser plotados diretamente após a geração, apenas executando o programa plot.py (escrito em Python)
   em um interpretador.   


### Passos para Executar o Projeto

1. Defina as variáveis presentes em 'main.f90';
2. Chame as funções com os parâmetros (h, tf, r0, phi0, vr0, vphi0, t, r, phi, vr, vphi);
3. Execute plot.py para ver o plot das órbitas.

