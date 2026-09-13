export const PRE_INP_KEY = "axia-pre-inp";
export const PRE_NAME_KEY = "axia-pre-name";

export function sendInpToSolver(inp: string, name: string) {
  sessionStorage.setItem(PRE_INP_KEY, inp);
  sessionStorage.setItem(PRE_NAME_KEY, name);
}

export function takeInpFromPreprocessor(): { inp: string; name: string } | null {
  const inp = sessionStorage.getItem(PRE_INP_KEY);
  if (!inp) return null;
  const name = sessionStorage.getItem(PRE_NAME_KEY) ?? "Präprozessor";
  sessionStorage.removeItem(PRE_INP_KEY);
  sessionStorage.removeItem(PRE_NAME_KEY);
  return { inp, name };
}
