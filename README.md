# Claude Code Skills

A small, personal collection of Codex/Claude skills. Each skill lives in its own folder with a `SKILL.md` that contains the trigger description, concise workflow, and any supporting references.

## Repository Layout
- `bioinformatics/` – domain‑specific skills (CLIP-seq, LACE-seq).
- `LICENSE` – MIT license for the written material in this repo.
- `.gitignore` – ignores local scratch files (e.g., `bioinformatics/pipeline.md` copies).

## Available Skills
- **clipseq-workflow** (`bioinformatics/clipseq-workflow/SKILL.md`): end‑to‑end CLIP/eCLIP processing from raw FASTQ to peaks/motifs—QC, adapter/UMI handling, alignment, UMI deduplication, peak calling (PEAKachu/PureCLIP), tracks, and motif/annotation tips.
- **laceseq-workflow** (`bioinformatics/laceseq-workflow/SKILL.md`): LACE/PARN-LACE-seq preprocessing and analysis—adapter/polyA trimming (cutadapt), Bowtie rRNA filtering and genome mapping, QC/coverage, repeat enrichment, metaplot generation with outlier filtering, and motif discovery.

## Using a Skill
1) Open the relevant `SKILL.md`; follow the Quick Start and core commands.
2) Load `references/` files only when needed (full pipelines, templates, SLURM examples).
3) For Codex/Claude runs, keep prompts tight—describe the task and name the skill to trigger it.
4) Validate changes with `scripts/quick_validate.py <skill>` (from the skill-creator toolkit).

## Contributing
- Add new skills in their own folder with a `SKILL.md`; include `agents/openai.yaml` if UI metadata is needed.
- Keep SKILL.md concise; move long command sets to `references/`.
- Prefer runnable snippets (bash/python) over prose when repeatable steps exist.
- Use ASCII unless the existing file requires otherwise.
- Run `quick_validate.py` before committing.
