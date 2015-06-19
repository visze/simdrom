package de.charite.compbio.simdrom.sampler;

import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.reference.IndexedFastaSequenceFile;
import htsjdk.variant.variantcontext.VariantContext;

import java.io.File;
import java.io.FileNotFoundException;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Iterator;
import java.util.List;
import java.util.Random;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import com.google.common.math.DoubleMath;

/**
 * @author <a href="mailto:max.schubach@charite.de">Max Schubach</a>
 *
 */
public class DeNovoSampler implements Iterator<VariantContext> {

	private double deNovoRate;
	private IndexedFastaSequenceFile referenceFile;
	private VariantContext nextVariant;
	ListMultimap<SAMSequenceRecord, Long> deNovoPositions;

	public DeNovoSampler(double deNovoRate, String referenceFile) throws FileNotFoundException {
		this.deNovoRate = deNovoRate;
		this.referenceFile = new IndexedFastaSequenceFile(new File(referenceFile));
		this.nextVariant = sampleNextVariant();
		calculateVariants();
	}

	private void calculateVariants() {
		deNovoPositions = ArrayListMultimap.<SAMSequenceRecord, Long> create();
		long size = referenceFile.getSequenceDictionary().getReferenceLength();
		Random random = new Random();
		List<Long> values = new ArrayList<Long>();
		for (double i = deNovoRate*(double) size; i > 0; i--) {
			if (i >=1)
				values.add(random.nextLong() % size);
			else {
				long value = random.nextLong() % DoubleMath.roundToLong(size * 1.0/i, RoundingMode.HALF_EVEN);
				if (value < size )
					values.add(value);
			}
		}
		Collections.sort(values);
		long position = 0;
		for (SAMSequenceRecord sequence : referenceFile.getSequenceDictionary().getSequences()) {
			int length = sequence.getSequenceLength();
			for (Long value : values) {
				if (value >= position && value < length) {
					deNovoPositions.put(sequence, value-position);
				}
				if (value >= length)
					break;
			}
			position += length;
		}
	}

	@Override
	public boolean hasNext() {
		return nextVariant != null;
	}

	@Override
	public VariantContext next() {
		VariantContext output = nextVariant;
		this.nextVariant = sampleNextVariant();
		return output;
	}

	private VariantContext sampleNextVariant() {
	}

}
