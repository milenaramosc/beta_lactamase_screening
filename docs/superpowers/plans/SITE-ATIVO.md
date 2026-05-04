Você é um agente de desenvolvimento especializado em Python, bioinformática estrutural e pipelines científicos em terminal.

Contexto:
O projeto possui um Target Profiler em:

src/betalactamase_engine/target_profiler.py
scripts/profile_target.py

O profiler analisa arquivos PDB de beta-lactamases, identifica SITE records, resíduos candidatos, metais, ligantes, águas e gera um target_profile.json.

Objetivo desta tarefa:
Adicionar uma etapa interativa e opcional para permitir que o usuário informe manualmente o sítio ativo quando o arquivo PDB não possuir registros SITE.

Importante:
Não implementar interface gráfica.
Não implementar sítios alostéricos.
Não implementar dinâmica molecular.
Não implementar docking novo.
Não implementar ADMET.
Não implementar algoritmo genético.

Requisito principal:
Quando o PDB não tiver SITE records, o script deve perguntar no terminal se o usuário sabe informar o sítio ativo.

Fluxo desejado:

1. O profiler deve primeiro tentar usar registros SITE do PDB.

2. Se existirem registros SITE:
   - usar normalmente o método atual;
   - detection_method deve continuar como "SITE_RECORD";
   - não perguntar nada ao usuário.

3. Se NÃO existirem registros SITE:
   - o script deve perguntar no terminal:

"O arquivo PDB não possui registros SITE anotados.
Você sabe informar os resíduos do sítio ativo? 
Digite no formato CHAIN:RESNAME:RESID separados por vírgula.
Exemplo: A:SER:70,A:LYS:73,A:GLU:166
Ou pressione ENTER para tentar detecção automática:"

4. Se o usuário digitar uma entrada válida:
   - usar esses resíduos como sítio ativo;
   - calcular o centro geométrico com base nos átomos desses resíduos;
   - definir:
     detection_method = "USER_PROVIDED_SITE";
   - adicionar esses resíduos em candidate_catalytic_residues;
   - adicionar em classification.evidence algo como:
     "Active site residues provided by user";
   - adicionar em quality_control.notes:
     "SITE records were absent; active site was provided manually by user."

5. Se o usuário pressionar ENTER sem informar nada:
   - não falhar;
   - seguir para a detecção automática atual.

6. Se o usuário digitar resíduos inválidos:
   - avisar no terminal quais resíduos não foram encontrados;
   - perguntar uma vez se deseja tentar novamente;
   - se continuar inválido ou vazio, seguir para detecção automática.

7. A detecção automática deve manter a hierarquia atual:
   - COCRYSTAL_LIGAND;
   - METAL_CENTER;
   - CATALYTIC_RESIDUES;
   - GEOMETRIC_FALLBACK.

8. Adicionar argumento opcional ao script:

--active-site

Esse argumento permite informar o sítio ativo diretamente sem prompt interativo.

Exemplo:

python scripts/profile_target.py \
  --pdb data/structures/betalac13.pdb \
  --chain A \
  --active-site A:SER:35,A:ASN:66,A:SER:86,A:SER:186 \
  --out results/experiments/teste/target_profile.json \
  --verbose

Se --active-site for fornecido:
   - usar esses resíduos;
   - não perguntar nada no terminal;
   - detection_method = "USER_PROVIDED_SITE".

9. Adicionar argumento opcional:

--no-interactive

Quando esse argumento for usado:
   - o script não deve perguntar nada ao usuário;
   - se não houver SITE, deve seguir direto para detecção automática.

Exemplo:

python scripts/profile_target.py \
  --pdb data/structures/betalac13.pdb \
  --chain A \
  --no-interactive \
  --out results/experiments/teste/target_profile.json \
  --verbose

10. O comportamento padrão deve ser:
   - se rodar no terminal normal e não houver SITE, perguntar ao usuário;
   - se --no-interactive for usado, não perguntar;
   - se --active-site for usado, usar o sítio informado.

11. Validar formato da entrada:

Formato aceito:
CHAIN:RESNAME:RESID

Exemplos válidos:
A:SER:70
A:LYS:73
A:GLU:166
A:SER:35,A:ASN:66,A:SER:86

Regras:
- CHAIN deve ser texto curto, geralmente A, B, C etc.
- RESNAME deve ter 3 letras.
- RESID deve ser número inteiro.
- aceitar espaços extras e remover automaticamente.
- converter RESNAME para maiúsculo.

12. Atualizar o target_profile.json.

Quando o usuário informar o sítio ativo, o JSON deve conter:

"active_site": {
  "center": [...],
  "detection_method": "USER_PROVIDED_SITE",
  "radius_angstrom": 8.0,
  "candidate_catalytic_residues": [...]
}

E também:

"quality_control": {
  "has_site_records": false,
  "used_fallback": false,
  "notes": [
    "SITE records were absent; active site was provided manually by user."
  ]
}

13. Se alguns resíduos informados não forem encontrados:
   - registrar em warnings;
   - usar apenas os resíduos válidos, se houver pelo menos um;
   - se nenhum for válido, seguir para detecção automática.

14. Atualizar docs/target_profiler.md explicando:

- o que acontece quando o PDB não tem SITE;
- como informar o sítio ativo interativamente;
- como usar --active-site;
- como usar --no-interactive;
- exemplos de comandos.

Critérios de aceite:

1. Com PDB que possui SITE:
   - usar SITE_RECORD;
   - não perguntar nada.

2. Com PDB sem SITE e sem --no-interactive:
   - perguntar ao usuário.

3. Com PDB sem SITE e usuário informando resíduos válidos:
   - usar USER_PROVIDED_SITE.

4. Com PDB sem SITE e usuário pressionando ENTER:
   - seguir para detecção automática atual.

5. Com --active-site:
   - usar diretamente os resíduos fornecidos;
   - não perguntar nada.

6. Com --no-interactive:
   - não perguntar nada;
   - seguir direto para detecção automática.

7. O pipeline antigo não pode quebrar.