// The page's Content-Security-Policy, for the tests: the <meta> in index.html and the header the service worker adds
// (sw.js) must be one policy, which lets in the page's one inline <style> by its hash.  policy() reads both from the
// source (csp.test.mjs) or from the site as published (site.mjs, after web/stamp.sh); watchViolations() records what
// the browser blocks while a test uses the page.
import { createHash } from 'node:crypto';

const directives = (policy) => new Map(policy.split(';').map((d) => d.trim().split(/\s+/)).filter((d) => d[0])
  .map(([name, ...sources]) => [name.toLowerCase(), sources]));

export function policy(page, sw) {
  const html = page.replace(/<!--[\s\S]*?-->/g, '');
  const metas = [...html.matchAll(/<meta http-equiv="Content-Security-Policy" content="([^"]*)">/gi)];
  if (metas.length !== 1) throw new Error(`index.html has ${metas.length} Content-Security-Policy <meta> tags, not one`);
  const meta = metas[0][1];
  const styles = [...html.matchAll(/<style>([\s\S]*?)<\/style>/g)];
  if (styles.length !== 1) throw new Error(`index.html has ${styles.length} <style> elements, not one`);
  // a policy in a <meta> covers only what comes after it
  const first = html.search(/<(style|script|link)\b/i);
  if (metas[0].index > first) throw new Error('the policy\'s <meta> comes after the first <style>, <script> or <link>');
  // a browser hashes the text as parsed, where a CRLF is read as LF
  const hash = `'sha256-${createHash('sha256').update(styles[0][1].replace(/\r\n?/g, '\n'), 'utf8').digest('base64')}'`;
  const styleSrc = directives(meta).get('style-src') ?? [];
  if (styleSrc.length !== 1 || styleSrc[0] !== hash) {
    throw new Error(`style-src is ${styleSrc.join(' ') || 'missing'}, but the <style> in index.html hashes to ${hash}`);
  }
  // inline scripts, event handler attributes and style attributes would be blocked
  const inline = [...html.matchAll(/<script\b(?![^>]*\ssrc=)[^>]*>|<[a-z][^>]*\s(?:on[a-z]+|style)=[^>]*>/gi)].map((m) => m[0]);
  if (inline.length) throw new Error(`index.html has markup the policy blocks: ${inline.join(' ')}`);
  // the service worker's header: the string literals that make up const POLICY
  const literal = /const POLICY = ((?:"[^"]*"\s*\+?\s*)+);/.exec(sw)?.[1];
  if (!literal) throw new Error('sw.js has no const POLICY = "…";');
  const header = [...literal.matchAll(/"([^"]*)"/g)].map((m) => m[1]).join('');
  if (header !== `${meta}; frame-ancestors 'none'`) {
    throw new Error(`the policy sw.js sends is not the page's plus frame-ancestors 'none':\n  sw.js:      ${header}\n  index.html: ${meta}`);
  }
  if (!/headers\.set\('Content-Security-Policy', POLICY\)/.test(sw)) throw new Error('sw.js does not send POLICY as a header');
  return { meta, header, hash };
}

// Every violation a page of `context` reports goes into `seen`: its documents' securitypolicyviolation events, the
// console messages about the policy, and in Chromium the log of the page's workers too (the engine's worker and its
// threads, whose console Playwright does not pass on).
const aboutPolicy = /content[- ]security[- ]policy/i;
export async function watchViolations(context, seen) {
  await context.exposeBinding('__sscViolation', (_, text) => seen.push(text));
  await context.addInitScript(() => {
    addEventListener('securitypolicyviolation', (e) => {
      self.__sscViolation?.(`${e.effectiveDirective} blocked ${e.blockedURI || '(inline)'} at ${e.sourceFile || e.documentURI}`
        + `:${e.lineNumber}${e.sample ? ` "${e.sample}"` : ''}`);
    });
  });
  context.on('console', (m) => { if (aboutPolicy.test(m.text())) seen.push(`console ${m.type()}: ${m.text()}`); });
  if (context.browser()?.browserType().name() !== 'chromium') return;
  context.on('page', async (page) => {
    try {
      const session = await context.newCDPSession(page);
      session.on('Log.entryAdded', ({ entry }) => {
        if (entry.source === 'worker' && aboutPolicy.test(entry.text)) seen.push(`worker log: ${entry.text}`);
      });
      await session.send('Log.enable');  // what was logged before is sent too
    } catch {}  // a page closed at once
  });
}
