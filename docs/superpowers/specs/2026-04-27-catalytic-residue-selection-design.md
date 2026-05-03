## Contexto e objetivo
O script `scripts/generate_inhibitors_from_protein.py` identifica residuos cataliticos com base em SITE records do PDB e, se nao existir, usa um array fixo `known_catalytic`. Esse fallback esta restritivo. O objetivo e permitir que o usuario escolha, em tempo de execucao, se deseja usar residuos conhecidos (informando quais) ou buscar em todo o PDB sem pre-conhecidos.

## Abordagens consideradas
1) Prompt interativo na CLI (recomendado): perguntar durante a execucao se deve usar residuos conhecidos e quais, ou buscar no PDB inteiro. Simples e alinhado ao uso atual.
2) Flags de CLI: `--catalytic-mode` e `--catalytic-residues`. Melhor para automatizacao, mas muda o uso do script.
3) Config em arquivo: mais flexivel, mas adiciona complexidade desnecessaria.

## Decisao
Usar prompt interativo para orientar a selecao de residuos cataliticos. Mantem o uso atual e remove a restricao do array fixo, sem adicionar novos parametros obrigatorios. Trade-offs: (1) reduz automatizacao direta em pipelines, (2) simplifica uso manual e reduz necessidade de lembrar flags. Em execucoes nao interativas (stdin nao-TTY ou EOF), o script assume `all` e avisa, evitando bloqueio. Se surgir necessidade de batch, uma evolucao futura pode adicionar flags opcionais sem quebrar o fluxo.

## Arquitetura e componentes
- `main()`: coleta a escolha do usuario e determina o fallback.
- `prompt_catalytic_selection()`: encapsula perguntas, validacao e retorno (`fallback_mode`, `fallback_resnames`).
- `analyze_binding_site()`: usa `identify_catalytic_residues()` com os parametros do fallback.
- `identify_catalytic_residues()`: aplica prioridade de SITE records e fallback escolhido.

## Design proposto
### Fluxo
- Antes de analisar o sitio ativo, o script pergunta:
  - Usar residuos cataliticos conhecidos? (s/n)
  - Se sim, solicita lista de residuos (ex.: `SER,LYS,GLU`), valida e usa essa lista.
  - Se nao, busca no PDB inteiro sem filtro por nomes de residuos.
- SITE records continuam tendo prioridade. Se existem, sao usados independentemente do modo escolhido (como hoje). O modo escolhido afeta apenas o fallback.

### Fluxo de dados
1) Entrada: `stdin` coleta escolha do modo e lista de residuos (quando aplicavel).
2) Validacao: normaliza para uppercase, filtra tokens vazios, valida tamanho 3 e apenas letras.
3) Decisao: define `fallback_mode` e `fallback_resnames`.
4) Selecoes: `identify_catalytic_residues` prioriza SITE records; se ausente, aplica fallback.
5) Saida: lista de atomos cataliticos para calcular centro e vizinhanca do sitio.

### Mudancas nas funcoes
- `identify_catalytic_residues(structure, site_residues, fallback_mode, fallback_resnames)`
  - `fallback_mode`: `"known"` ou `"all"`.
  - `fallback_resnames`: lista de nomes (ex.: `["SER","LYS"]`), usada apenas em `known`.
  - Se `site_residues` existir, usa-os e retorna.
  - Se `fallback_mode == "known"`, filtra residuos por `fallback_resnames`.
  - Se `fallback_mode == "all"`, inclui atomos de todos os residuos do modelo (apenas aminoacidos padrao, excluindo agua/HOH/WAT e ligantes hetero).
  - Se `site_residues` existir mas estiver vazio/invalidado, aplica fallback normalmente.

### Validacao
- Lista de residuos: uppercase, 3 letras, separados por virgula.
- Se `known` e lista vazia, pedir novamente.
- Rejeitar tokens com caracteres nao alfabeticos ou tamanho diferente de 3.
- Aceitar entradas `s`, `sim`, `y`, `yes` como afirmativo; `n`, `nao`, `no` como negativo.
 - Para entradas desconhecidas (ex.: `talvez`), repetir prompt ate 3 vezes.

### Erros e mensagens
- Mensagem clara quando `site_residues` nao existe e o fallback e ativado.
- Mensagens com o modo escolhido (known com lista ou all).
- Em stdin nao-TTY, EOF ou interrupt durante prompt: usar fallback padrao `all` e imprimir aviso.
- Em entrada invalida: repetir prompt (max 3 tentativas), depois usar `all` e avisar.

## Testes manuais
1) PDB com SITE records: garantir que SITE e usado e o fallback escolhido nao altera o resultado.
2) PDB sem SITE records:
   - Modo `known` com lista pequena (ex.: `SER,LYS`) e confirmar que residuos identificados refletem isso.
   - Modo `all` e confirmar que o site e calculado usando todos os residuos (apenas aminoacidos padrao).
3) Entrada invalida:
   - Lista vazia ou tokens invalidos e confirmar retry e fallback para `all` apos 3 tentativas.
4) Execucao nao interativa (EOF): confirmar aviso e uso de `all`.
5) Regressao:
   - Confirmar que SITE records, quando presentes e validos, sempre prevalecem.

## Integracao no script
- Inserir prompt no `main()` antes de `analyze_binding_site`.
- Passar `fallback_mode` e `fallback_resnames` para `analyze_binding_site`, que repassa a `identify_catalytic_residues`.
- Manter comportamento atual quando SITE records existem.
- Considerar multiplos modelos: manter uso do `model = structure[0]` como hoje.
- Tornar o prompt opcional: se stdin nao for TTY, pular prompt e usar `all`.
