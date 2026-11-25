"""
Exact count for drafts/route_twin_loop_by_hand.md: the routes of III.15a (A) (walk read in t, tick = shift then coin
R = cos th 1 + sin th T), N ticks from S = (0,0) in the forward sense, on an unbounded line. Weights are products of the
signed coin entries; nothing is solved or approximated.
"""
import itertools, sys
import numpy as np
ok_all = True
def ok(msg, cond):
    global ok_all; ok_all &= bool(cond); print(("PASS  " if cond else "FAIL  ") + msg); sys.stdout.flush()
th = 0.4; c, s = np.cos(th), np.sin(th)
def coin(new, old):                             # R[new, old] with senses +1 -> row 0, -1 -> row 1
    R = np.array([[c, s], [-s, c]]); return R[0 if new > 0 else 1, 0 if old > 0 else 1]
def route(senses, launch=+1):                   # senses carried on each tick; the coin after each tick sets the next; the last sets the arrival sense
    x, t, w, pts, old = 0, 0, 1.0, [(0, 0)], launch
    for e in senses:
        x += e; t += 1; pts.append((x, t)); w *= coin(e, old); old = e
    return pts, w
def all_routes(N, launch=+1):
    out = {}
    for sn in itertools.product((+1, -1), repeat=N):
        pts, w = route(sn, launch); out.setdefault((pts[-1], sn[-1]), []).append((sn, w, pts))
    return out
def back_weight(senses, arrival):               # III.15a (B): the return along a route reversed = the same coin entries, other order
    s_ = list(senses) + [arrival]
    return float(np.prod([coin(s_[i + 1], s_[i]) for i in range(len(senses))]))

for d1, d2 in ((2, 1), (3, 2)):
    N = d1 + d2
    r = (+1,)*d1 + (-1,)*d2                      # the route: forward d1, then backward d2; arrival sense -1
    tw = (-1,)*d2 + (+1,)*d1                     # its twin: backward d2, then forward d1; arrival sense +1
    pts_r, w_r = route(r, +1); pts_tw, w_tw_raw = route(tw, -1); w_tw = -w_tw_raw     # launched with the sign of T: v(0) = -e_- g
    Q = pts_r[-1]
    # the twin reversed as a return from Q: its own entries in the other order (launch sense -1 at S, arrival +1 at Q)
    w_back = back_weight(tw, +1)
    pair = w_r*w_back
    ok("d1 = %d, d2 = %d: route and twin meet at Q = %s (twin ends at %s); w_route = %+.5f = -sin cos^%d = %+.5f; "
       "w_twin (second column, launched with -1) = %+.5f = w_route; return along the twin = %+.5f = +sin cos^%d; "
       "route out x twin back = %+.5f = -sin^2 cos^%d = %+.5f"
       % (d1, d2, Q, pts_tw[-1], w_r, N - 1, -s*c**(N - 1), w_tw, w_back, N - 1, pair, 2*N - 2, -s*s*c**(2*N - 2)),
       pts_tw[-1] == Q and abs(w_r + s*c**(N - 1)) < 1e-15 and abs(w_tw - w_r) < 1e-15 and abs(w_back - s*c**(N - 1)) < 1e-15
       and abs(pair + s*s*c**(2*N - 2)) < 1e-15)
    # the two families of closed figures at S: out from (S,+) with every return; returns to (S,+) are reversed out routes
    # launched +, returns to (S,-) are reversed out routes launched -, each with its own out weight (III.15a (B), (C)2)
    Rp = all_routes(N, +1); Rm = all_routes(N, -1)
    same = sum(sum(w for _, w, _ in v)**2 for v in Rp.values())
    other = sum(sum(w for _, w, _ in Rp[k])*sum(w for _, w, _ in Rm.get(k, [])) for k in Rp)
    Qm = sum(w for _, w, _ in Rp.get((Q, -1), [])); ev = sum(wi*wj for _, wi, _ in Rp.get((Q, -1), []) for _, wj, _ in Rp.get((Q, -1), []))
    ok("   route arrives at Q in sense %+d, twin in sense %+d (same point, opposite counters); event weight at (Q,-) for source (S,+) = %.6f = Q_-(Q)^2 "
       "(III.15a (C)3, the twin reversed is not among its returns); the family out from (S,+) home to (S,+) over every far point = %.12f (VI.1 (A)); "
       "the family out from (S,+) home to (S,-) = %.1e: the route-twin figure (%+.5f) is one term of the family that totals zero"
       % (r[-1], tw[-1], ev, same, other, pair), r[-1] == -1 and tw[-1] == +1 and abs(ev - Qm*Qm) < 1e-12 and abs(same - 1) < 1e-12 and abs(other) < 1e-12)
    # (added 10 Nov 2026) the text's "route x twin at the meeting point is -v_+ u_-, a term of W = u_+ v_- - v_+ u_-":
    # u = the column of all routes launched forward, v = the column of all routes launched backward with the sign of T
    u_minus = Rp.get((Q, -1), []); v_plus = Rm.get((Q, +1), [])
    in_u = any(sn == r and abs(w - w_r) < 1e-15 for sn, w, _ in u_minus); in_v = any(sn == tw and abs(-w - w_tw) < 1e-15 for sn, w, _ in v_plus)
    U_p = sum(w for _, w, _ in Rp.get((Q, +1), [])); U_m = sum(w for _, w, _ in u_minus)
    V_p = -sum(w for _, w, _ in v_plus); V_m = -sum(w for _, w, _ in Rm.get((Q, -1), []))
    Wq = U_p*V_m - V_p*U_m
    cross = sum((-wv)*wu for _, wv, _ in v_plus for _, wu, _ in u_minus)          # all the terms of v_+ u_- at Q
    ok("   the route is one of the routes making u_-(Q) and the twin one of those making v_+(Q); their product %+.5f is one of the %d terms of v_+ u_- "
       "(= %+.6f), which enters W = u_+ v_- - v_+ u_- = %+.6f with a minus: the figure is a term of what turns the phase"
       % (w_tw*w_r, len(v_plus)*len(u_minus), cross, Wq), in_u and in_v and abs(cross - V_p*U_m) < 1e-12 and abs(w_tw*w_r - w_r**2) < 1e-15)
    # VI.2 (C): uniform field, angle E x per tick at the site reached; route minus twin = E x enclosed area 2 d1 d2
    E = 0.37
    ang_r = sum(E*x for x, _ in pts_r[1:]); ang_tw = sum(E*x for x, _ in pts_tw[1:])
    ok("   uniform field, angle E x per tick: route %.4f - twin %.4f = %.4f = E x (area 2 d1 d2) = %.4f (VI.2 (C): the figure reads the enclosed flux)"
       % (ang_r, ang_tw, ang_r - ang_tw, E*2*d1*d2), abs(ang_r - ang_tw - E*2*d1*d2) < 1e-12)
print("ALL PASS" if ok_all else "FAIL"); raise SystemExit(0 if ok_all else 1)
