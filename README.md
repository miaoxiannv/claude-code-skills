# Claude Code Skills

A small, personal collection of Codex/Claude skills. Each skill lives in its own folder with a `SKILL.md` that contains the trigger description, concise workflow, and any supporting references.

## Repository Layout
- `bioinformatics/` – domain‑specific skills; currently contains a CLIP‑seq workflow guide.
- `LICENSE` – MIT license for the written material in this repo.

## Available Skills
- **clipseq-workflow** (`bioinformatics/clipseq-workflow/SKILL.md`): end‑to‑end CLIP/eCLIP processing from raw FASTQ to peaks/motifs—QC, adapter/UMI handling, alignment, UMI deduplication, peak calling (PEAKachu/PureCLIP), tracks, and motif/annotation tips.

## Using a Skill
1) Open the relevant `SKILL.md` to follow the outlined steps and command templates.
2) If the skill references files under `references/`, load only what you need (e.g., pipeline command snippets).
3) Keep edits minimal and scoped to the skill you are using to reduce context bloat.

## Contributing
- Add new skills in their own folder with a `SKILL.md` and any minimal `references/` assets.
- Prefer concise checklists over long narratives; link out for deep background when needed.
- Stay within ASCII unless the existing file already uses other characters.
