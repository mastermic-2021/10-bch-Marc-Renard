/* x is a generator of Fq obtained via ffgen(q) */
[q,t,m] = readvec("input.txt");
f = ffgen(q,'a); \\création de f, générateur de F128
encodefq(i,x)=subst(Pol(digits(i,x.p)),'x,x);
decodefq(z,x)=fromdigits(Vec(z.pol),x.p);
int2fqx(i,x)=Polrev([encodefq(z,x)|z<-digits(i,x.p^x.f)]);
fqx2int(p,x)=fromdigits([decodefq(z,x)|z<-Vecrev(p)],x.p^x.f);



\\ Le message a été caché en prenant un mot de code valide d'un code BCH, et en y ajoutant le message ce qui a créé des erreurs dans le code
\\On cherche donc à extraire cette erreur ajoutée pour avoir le message.
\\Nous allons pour cela calculer le syndrôme, et utiliser les approximants de Padé pour trouver:
\\	- Le polynôme E localisateur de l'erreur
\\	- Le polynôme R qui lui, nous servira à calculer la valeur de ces erreurs (E intervient aussi dans ce calcul) 

y=int2fqx(m, f); \\transforme le message qui est un nombre en un polynome de Fq


syndrome(y, a, c) = sum(i=0, 2 * t - 1, subst(y, 'x, a^(i+c)) * 'x^i); \\définition du syndrôme tel que définit dans le cours
\\On sait que les racines sont 2t-1 puissances de a consécutives, mais on ne sait pas à partir de quelles puissance ces puissances débutent, d'où la constante c.
reveleMessage() = {
	\\Boucle de taille 128, je suppose que le retour de ffprimeroot suit un schéma cyclique et comme le cardinal de F128 est 128 on va parcourrir toutes les valeurs possibles de retour de ffprimroot
	for(k=0, 128,
		a = ffprimroot(f); \\ a est donc une racine 128ème primitive de l'unité. Cette ligne est mise dans la boucle car elle ne donne pas toujours le même résultat
		for(j=0, 127,
			ratio=bestapprPade(Mod(syndrome(y, a, j), x^(2*t)));
			R = numerator(ratio);
			E = denominator(ratio); \\Polynôme localisateur de l'erreur
			erreur = List(); \\liste qui va contenir la liste des erreurs commises
			for(i = 0, 126, pol = subst(E, 'x, a^(-i)); if (pol == 0, val = subst(R/deriv(E) * (x^(j-1)), 'x, a^(-i)); l = fqx2int(val, f); listput(erreur, l)));
			\\le message a été caché en insérant les erreurs, on peut donc supposer que l'on doit trouver un nombre d'erreur assez grand pour représenter un message
			if(#erreur > 10 ,print(Strchr(Vec(erreur))); return);
		);
	);
};


reveleMessage();
