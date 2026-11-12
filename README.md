# Tucanos documentation

This repository contains the mdBook sources for the `tucanos` project:

- Project: [tucanos](https://github.com/tucanos/tucanos)
- Published docs: [https://tucanos.github.io](https://tucanos.github.io)
- Documentation engine: [mdBook](https://rust-lang.github.io/mdBook/)

## Quick start

Build the book locally:

```bash
mdbook build
```

Serve locally with live reload:

```bash
mdbook serve --open
```

## Testing documentation examples

### Rust snippets

Rust code fences are validated with `mdbook test`.

```bash
mdbook test .
```

When snippets import `tmesh` or `tucanos`, make sure those crates are available in library search paths. A common setup is:

```bash
mdbook test . -L <path>/target/debug/deps -L <path>/target/debug
```

### Python snippets

Python examples are currently documentation examples only and are not executed automatically by mdBook CI.