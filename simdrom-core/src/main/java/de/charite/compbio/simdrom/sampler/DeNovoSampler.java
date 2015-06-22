package de.charite.compbio.simdrom.sampler;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.samtools.reference.ReferenceSequence;
import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;

import java.io.File;
import java.io.FileNotFoundException;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.Random;

import com.google.common.collect.ImmutableSet;
import com.google.common.math.DoubleMath;
import com.google.common.primitives.Ints;

/**
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class DeNovoSampler implements Iterator<VariantContext> {

	private double deNovoRate;
	private IndexedFastaSequenceFile referenceFile;
	Queue<VariantContext> deNovoPositions;

	public DeNovoSampler(double deNovoRate, String referenceFile) throws FileNotFoundException {
		this.deNovoRate = deNovoRate;
		this.referenceFile = new IndexedFastaSequenceFile(new File(referenceFile));
		calculateVariants();
	}

	private void calculateVariants() {

		// sample positions
		long size = referenceFile.getSequenceDictionary().getReferenceLength();
		Random random = new Random();
		List<Long> values = new ArrayList<Long>();
		for (double i = deNovoRate * (double) size; i > 0; i--) {
			if (i >= 1)
				values.add(Math.abs(random.nextLong()) % size);
			else {
				long value = Math.abs(random.nextLong()) % DoubleMath.roundToLong(size * 1.0 / i, RoundingMode.HALF_EVEN);
				if (value < size)
					values.add(value);
			}
		}

		// get variant
		Collections.sort(values);
		long position = 0;
		for (SAMSequenceRecord sequence : referenceFile.getSequenceDictionary().getSequences()) {
			int length = sequence.getSequenceLength();
			for (Long value : values) {
				if (value >= position && value < length+position) {
					String chr = sequence.getSequenceName();
					int pos = Ints.checkedCast(value - position);
					ReferenceSequence refSeq = referenceFile.getSubsequenceAt(chr, pos, pos);
					Allele ref = Allele.create(refSeq.getBases(), true);
					Allele alt = createNewAllele(Allele.create(refSeq.getBases()));
					VariantContext vc = new VariantContextBuilder("deNovo", chr, value - position, value - position,
							ImmutableSet.<Allele> of(ref, alt)).make();
					getDeNovoPositions().add(vc);
				}
				if (value >= length+position)
					break;
			}
			position += length;
		}
	}

	private Allele createNewAllele(Allele b) {
		Random random = new Random();
		int i = random.nextInt(3);
		Allele n = getNucleotide(i);
		if (n == b)
			n = getNucleotide(3);
		return n;
	}

	private Allele getNucleotide(int i) {
		switch (i) {
		case 0:
			return Allele.create(new byte[]{'A'});
		case 1:
			return Allele.create(new byte[]{'T'});
		case 2:
			return Allele.create(new byte[]{'C'});
		default:
			return Allele.create(new byte[]{'G'});
		}
	}

	@Override
	public boolean hasNext() {
		return getDeNovoPositions().peek() != null;
	}

	@Override
	public VariantContext next() {
		return getDeNovoPositions().poll();
	}

	private Queue<VariantContext> getDeNovoPositions() {
		if (deNovoPositions == null)
			deNovoPositions = new LinkedList<VariantContext>();
		return deNovoPositions;
	}
	@Override
	public void remove() {
		throw new UnsupportedOperationException();
	}

}
