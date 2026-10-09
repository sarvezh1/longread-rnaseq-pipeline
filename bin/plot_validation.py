#!/usr/bin/env python3
"""Create a deterministic offline SVG line/scatter plot from a benchmark TSV."""

import argparse, csv, html
from collections import defaultdict
from pathlib import Path

COLORS=["#1b6ca8","#c44e52","#4c956c","#7a5195","#d17c0f","#4c4c4c"]


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--input",type=Path,required=True); parser.add_argument("--x",required=True); parser.add_argument("--y",required=True); parser.add_argument("--group",required=True); parser.add_argument("--title",required=True); parser.add_argument("--x-label",required=True); parser.add_argument("--y-label",required=True); parser.add_argument("--output",type=Path,required=True); args=parser.parse_args()
    with args.input.open(newline="",encoding="utf-8") as handle: source=list(csv.DictReader(handle,delimiter="\t"))
    points=defaultdict(list)
    for row in source:
        try: points[row[args.group]].append((float(row[args.x]),float(row[args.y])))
        except (KeyError,ValueError): continue
    values=[point for group in points.values() for point in group]
    if not values: raise SystemExit("No finite numeric points available")
    xmin,xmax=min(x for x,_ in values),max(x for x,_ in values); ymin,ymax=min(y for _,y in values),max(y for _,y in values)
    if xmin==xmax: xmin-=.5; xmax+=.5
    if ymin==ymax: ymin-=.5; ymax+=.5
    left,top,width,height=85,55,700,430
    sx=lambda x:left+(x-xmin)/(xmax-xmin)*width; sy=lambda y:top+height-(y-ymin)/(ymax-ymin)*height
    shapes=[]; legend=[]
    for index,(group,series) in enumerate(sorted(points.items())):
        color=COLORS[index%len(COLORS)]; ordered=sorted(series); coordinates=" ".join(f"{sx(x):.2f},{sy(y):.2f}" for x,y in ordered)
        shapes.append(f'<polyline points="{coordinates}" fill="none" stroke="{color}" stroke-width="2"/>'+"".join(f'<circle cx="{sx(x):.2f}" cy="{sy(y):.2f}" r="4" fill="{color}"/>' for x,y in ordered))
        legend.append(f'<rect x="810" y="{70+index*25}" width="14" height="14" fill="{color}"/><text x="832" y="{82+index*25}" font-size="13">{html.escape(group)}</text>')
    ticks=[]
    for index in range(6):
        fraction=index/5; x=xmin+fraction*(xmax-xmin); y=ymin+fraction*(ymax-ymin)
        ticks.append(f'<line x1="{sx(x):.2f}" y1="{top+height}" x2="{sx(x):.2f}" y2="{top+height+6}" stroke="#333"/><text x="{sx(x):.2f}" y="{top+height+23}" text-anchor="middle" font-size="12">{x:.3g}</text><line x1="{left-6}" y1="{sy(y):.2f}" x2="{left}" y2="{sy(y):.2f}" stroke="#333"/><text x="{left-10}" y="{sy(y)+4:.2f}" text-anchor="end" font-size="12">{y:.3g}</text>')
    svg=f'''<svg xmlns="http://www.w3.org/2000/svg" width="1050" height="570" viewBox="0 0 1050 570"><rect width="100%" height="100%" fill="white"/><text x="525" y="28" text-anchor="middle" font-size="20" font-weight="bold">{html.escape(args.title)}</text><line x1="{left}" y1="{top}" x2="{left}" y2="{top+height}" stroke="#333"/><line x1="{left}" y1="{top+height}" x2="{left+width}" y2="{top+height}" stroke="#333"/>{''.join(ticks)}{''.join(shapes)}{''.join(legend)}<text x="{left+width/2}" y="550" text-anchor="middle" font-size="15">{html.escape(args.x_label)}</text><text x="20" y="{top+height/2}" text-anchor="middle" font-size="15" transform="rotate(-90 20 {top+height/2})">{html.escape(args.y_label)}</text></svg>\n'''
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(svg,encoding="utf-8")


if __name__=="__main__": main()
